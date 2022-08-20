-- | Infer network de novo from data
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network.DeNovo
    ( saveAssociations
    , mkAssociations
    , mkNetwork
    , readNetwork
    , readAssociations
    , getTFBS
    , outputBindingEdges
    , outputCombinedEdges
    , outputNodes
    ) where

import Control.Arrow ((&&&))
import Control.Monad.State.Strict
import           Bio.Data.Bed hiding (NarrowPeak)
import           Bio.Data.Bed (NarrowPeak)
import Data.Conduit.Internal (zipSinks)
import qualified Data.IntervalMap.Strict as IM
import qualified Data.Set as S
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import qualified Data.HashMap.Strict                   as M
import qualified Data.Text                         as T
import IGraph

import Taiji.Utils
import Taiji.Core.RegulatoryElement
import           Taiji.Prelude

-- | Construct and save nodes and edges.
saveAssociations :: ( ATACSeq S (File '[Gzip] 'Bed)
                    , Either (File '[] 'NarrowPeak) (File '[Gzip] 'NarrowPeak)
                    , Maybe (File '[ChromosomeLoop] 'Bed)
                    , Maybe (File '[] 'Tsv) )  -- ^ (TFBS, promoter activity, HiC loops, Expression)
                 -> ReaderT TaijiConfig IO
                        (ATACSeq S (File '[] 'Other, File '[] 'Other))
saveAssociations (tfFl, peakFl, hicFl, expr) = do
    dir <- asks
        ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir)
        >>= getPath
    anno <- fromJust <$> asks _taiji_annotation
    let netEdges = dir ++ "/edges_combined.csv"
        netNodes = dir ++ "/nodes.csv"
        bindingEdges = dir ++ "/edges_binding.csv"
    liftIO $ do
        openSites <- case peakFl of
            Left fl -> readBed $ fl ^.location
            Right fl -> runResourceT $ runConduit $
                streamBedGzip (fl^.location) .| sinkList
        expr' <- (fmap . fmap) (\(a,b) -> (sqrt a, exp b)) $ case expr of
            Nothing -> return M.empty
            Just e -> readExpression 1 (B.pack $ T.unpack grp ) $ e^.location
        tfbs <- runResourceT $ runConduit $
            streamBedGzip (tfFl^.replicates._2.files.location) .|
            getTFBS (mkPeakMap openSites)
        promoters <- findActivePromoters openSites <$> readPromoters anno 

        let proc = loops .| findTargets tfbs promoters .|
                mkAssociations expr' .|
                zipSinks (outputCombinedEdges netEdges)
                (outputBindingEdges bindingEdges)
            loops = case hicFl of
                Nothing -> return ()
                Just fl -> read3DContact $ fl^.location
        runResourceT (execStateT (runConduit proc) S.empty) >>=
            outputNodes netNodes
    return $ tfFl & replicates.mapped.files .~
        ( emptyFile & location .~ netNodes
        , emptyFile & location .~ netEdges )
  where
    grp = tfFl^.groupName._Just
{-# INLINE saveAssociations #-}

getTFBS :: Monad m
        => BEDTree PeakAffinity  -- ^ Potential regulatory regions and its affinity scores
        -> ConduitT BED o m (BEDTree [SiteInfo])
getTFBS peaks = concatMapC f .| sinkList >>=
    return . (fmap . fmap) nub' . bedToTree (++)
  where
    f site = case IM.elems (within peaks site) of
        [] -> Nothing
        xs -> Just $ (bed, [SiteInfo (getTFName site) siteSc (maximum xs)])
      where
        bed = asBed (site^.chrom) (site^.chromStart) (site^.chromEnd) :: BED3
        siteSc = toSiteAffinity (fromJust $ site^.score)
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    nub' = M.elems . M.fromListWith (\a b -> if g a > g b then a else b) .
        map (\x -> (_tf_name x, x))
      where
        g = getSiteAffinity . _site_affinity

-- | Construct nodes and edges.
mkAssociations :: Monad m
               => M.HashMap GeneName (Double, Double)   -- ^ edge weight and node weight
               -> ConduitT (GeneName, ([TFBS], [TFBS]))
                           NetEdge
                           (StateT (S.Set NetNode) m) ()
mkAssociations expr = concatMapMC $ \(geneName, (ps, es)) -> do
    let edgeEnhancer = mkEdges geneName "enhancer" es
        edgePromoter = mkEdges geneName "promoter" ps
        (geneExpr, scaledGeneExpr) = M.lookupDefault (0.1, 1) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_weight = scaledGeneExpr
                           , _node_expression = Just geneExpr }
    modify' $ S.insert geneNode
    edgeCombined <- forM (groupEdgeByTF $ edgeEnhancer ++ edgePromoter) $ \xs -> do
        let (tfExpr, scaledTfExpr) = M.lookupDefault (0.1, 1) tfName expr
            tfNode = NetNode { _node_name = tfName
                             , _node_weight = scaledTfExpr
                             , _node_expression = Just tfExpr }
            tfName = _edge_from $ head xs
            wCombined = lp 2 $ map (_edge_binding_affinity . _edge_type) xs
        modify' $ S.insert tfNode
        return $ NetEdge { _edge_from = tfName
                         , _edge_to = geneName
                         , _edge_type = Combined (wCombined * tfExpr) }
    return $ edgePromoter ++ edgeEnhancer ++ edgeCombined
  where
    mkEdges geneName anno = filter
            ((>=edge_weight_cutoff) . _edge_binding_affinity . _edge_type) .
            map siteToEdge
      where
        siteToEdge site = NetEdge
            { _edge_from = _tf_name $ site^._data
            , _edge_to = geneName
            , _edge_type = Binding
                { _edge_binding_locus = convert site
                , _edge_binding_annotation = anno
                , _edge_binding_affinity = getEdgeWeight site } }
    groupEdgeByTF = groupBy ((==) `on` _edge_from) . sortBy (comparing _edge_from)
    getEdgeWeight x = sqrt $ siteSc * peakSc
      where
        siteSc = getSiteAffinity $ _site_affinity $ x^._data
        peakSc = getPeakAffinity $ _peak_affinity $ x^._data
{-# INLINE mkAssociations #-}

mkNetwork :: Monad m
          => M.HashMap GeneName (Double, Double)   -- ^ Edge weight and node weight
          -> ConduitT (GeneName, ([TFBS], [TFBS])) o m (Graph 'D NetNode Double)
mkNetwork expr = fmap fromLabeledEdges $ concatMapC mkEdges .| sinkList
  where
    mkEdges (geneName, (ps, es)) = flip map tfGroup $ \tfs ->
        let (edgeW, tfWeight) = M.lookupDefault (0.1, 1) tfName expr
            tfNode = NetNode { _node_name = tfName
                            , _node_weight = tfWeight
                            , _node_expression = Just edgeW }
            tfName = fst $ head tfs
            wCombined = lp 2 $ map snd tfs
        in ((geneNode, tfNode), wCombined * edgeW )
      where
        (geneExpr, geneWeight) = M.lookupDefault (0.1, 1) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_weight = geneWeight
                           , _node_expression = Just geneExpr }
        tfGroup = groupBy ((==) `on` fst) $ sortBy (comparing fst) $
            filter ((>=edge_weight_cutoff) . snd) $ map f $ ps ++ es
          where
            f site = (_tf_name $ site^._data, getEdgeWeight site)
            getEdgeWeight x = sqrt $ siteSc * peakSc
              where
                siteSc = getSiteAffinity $ _site_affinity $ x^._data
                peakSc = getPeakAffinity $ _peak_affinity $ x^._data
{-# INLINE mkNetwork #-}


--------------------------------------------------------------------------------
-- IO related functions
--------------------------------------------------------------------------------

-- | Save the edge information to files.
outputBindingEdges :: MonadResource m
                   => FilePath -> ConduitT NetEdge Void m ()
outputBindingEdges output = filterC isBinding .|
    (yield header >> mapC edgeToLine) .| unlinesAsciiC .| sinkFile output
  where
    header = ":START_ID,:END_ID,chr,start:int,end:int," <>
        "annotation,affinity,:TYPE"
    isBinding e = case _edge_type e of
        Binding{} -> True
        _ -> False
{-# INLINE outputBindingEdges #-}

-- | Save the edge information to files.
outputCombinedEdges :: MonadResource m
                    => FilePath -> ConduitT NetEdge Void m ()
outputCombinedEdges output = filterC isCombined .|
    (yield header >> mapC edgeToLine) .| unlinesAsciiC .| sinkFile output
  where
    header = ":START_ID,:END_ID,weight,:TYPE"
    isCombined e = case _edge_type e of
        Combined _ -> True
        _ -> False
{-# INLINE outputCombinedEdges #-}

-- | Save the node information to a file.
outputNodes :: FilePath -> S.Set NetNode -> IO ()
outputNodes fl = B.writeFile fl . B.unlines . (nodeHeader:) .
    map nodeToLine . S.toList
  where
    nodeHeader = "geneName:ID,expression,expressionZScore"
{-# INLINE outputNodes #-}

-- | Build the network from files containing the information of nodes and edges.
readNetwork :: FilePath  -- ^ nodes
            -> FilePath  -- ^ edges
            -> IO (Graph 'D NetNode Double)
readNetwork nodeFl edgeFl = do
    nodeMap <- M.fromList . map ((_node_name &&& id) . nodeFromLine) .
        tail . B.lines <$> B.readFile nodeFl
    runResourceT $ fromLabeledEdges' edgeFl (toEdge nodeMap)
  where
    toEdge nodeMap fl = sourceFileBS fl .| linesUnboundedAsciiC .|
        (dropC 1 >> mapC f) 
      where
        f l = case B.split ',' l of
            [f1,f2,f3,_] ->
              ( ( M.lookupDefault undefined (mk f2) nodeMap
                , M.lookupDefault undefined (mk f1) nodeMap )
                , readDouble f3 )
            _ -> error $ "Unexpected line: " <> show l
{-# INLINE readNetwork #-}

-- | Read network files as nodes and edges
readAssociations :: FilePath   -- ^ nodes
                 -> FilePath   -- ^ edges
                 -> IO ([NetNode], [((GeneName, GeneName), Double)])
readAssociations nodeFl edgeFl = do
    nds <- map nodeFromLine . tail . B.lines <$> B.readFile nodeFl
    es <- map f . tail . B.lines <$> B.readFile edgeFl
    return (nds, es)
  where
    f l = ( ( mk f2, mk f1), readDouble f3 )
        where
        [f1,f2,f3,_] = B.split ',' l
{-# INLINE readAssociations #-}

-- | Construct peak map from narrowpeaks.
mkPeakMap :: [NarrowPeak] -> BEDTree PeakAffinity
mkPeakMap = bedToTree max . map f
  where
    f x = ( asBed (x^.chrom) (center - 50) (center + 50) :: BED3
          , toPeakAffinity $ fromMaybe 5 $ x^.npPvalue )
      where
        center = case x^.npPeak of
          Nothing -> (x^.chromStart + x^.chromEnd) `div` 2
          Just c -> x^.chromStart + c
{-# INLINE mkPeakMap #-}