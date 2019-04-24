-- | Infer network de novo from data

{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network.DeNovo
    ( createLinkage
    , mkNetwork
    , readNodesAndEdges
    ) where

import Control.Arrow ((&&&))
import           Bio.Utils.Misc                    (readDouble)
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (getPath, asDir)
import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Conduit
import Data.Conduit.Internal (zipSinks)
import qualified Data.IntervalMap.Strict as IM
import           Control.Lens                      hiding (pre, to)
import qualified Data.Set as S
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import Data.Ord (comparing)
import           Data.List                         (foldl', maximumBy)
import qualified Data.HashMap.Strict                   as M
import           Data.Maybe                        (fromJust, mapMaybe)
import qualified Data.Text                         as T
import IGraph
import           Scientific.Workflow               hiding (_data)

import Taiji.Core.Network.Utils
import Taiji.Core.RegulatoryElement (findTargets)
import           Taiji.Types
import           Taiji.Constants (edge_weight_cutoff)

createLinkage :: ( ATACSeq S ( File tag1 'Bed         -- ^ Active promoters
                             , File tag2 'Bed         -- ^ TFBS
                             )
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ promoter activity
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ enhancer activity
                 , Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
                 , Maybe (File '[] 'Tsv)          -- ^ Expression
                 )
              -> WorkflowConfig TaijiConfig
                    (ATACSeq S (File '[] 'Other, File '[] 'Other))
createLinkage (atac, pro, enh, hic, expr) = do
    dir <- asks
        ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir)
        >>= getPath
    let netEdges = dir ++ "/edges_combined.csv"
        netNodes = dir ++ "/nodes.csv"
        bindingEdges = dir ++ "/edges_binding.csv"
    liftIO $ do
        expr' <- case expr of
            Nothing -> return M.empty
            Just e -> readExpression 1 (B.pack $ T.unpack grp ) $ e^.location
        activityPro <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_pro
        activityEnh <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_enh
        let proc = findTargets active_pro tfbs hic .|
                createLinkage_ activityPro activityEnh expr' .|
                sinkLinkage netEdges bindingEdges
        runResourceT (execStateT (runConduit proc) S.empty) >>=
            outputNodes netNodes
    return $ atac & replicates.mapped.files .~
        ( emptyFile & location .~ netNodes
        , emptyFile & location .~ netEdges )
  where
    fl_pro = either (^.location) (^.location) pro
    fl_enh = either (^.location) (^.location) enh
    (active_pro, tfbs) = atac^.replicates._2.files
    grp = atac^.groupName._Just
{-# INLINE createLinkage #-}

outputNodes :: FilePath -> S.Set NetNode -> IO ()
outputNodes fl = B.writeFile fl . B.unlines . (nodeHeader:) .
    map nodeToLine . S.toList
  where
    nodeHeader = "geneName:ID,expression,expressionZScore"
{-# INLINE outputNodes #-}

sinkLinkage :: MonadResource m
            => FilePath
            -> FilePath
            -> ConduitT NetEdge Void (StateT (S.Set NetNode) m) ()
sinkLinkage combined binding = zipSinks outputCombined outputBindSites >> return ()
  where
    outputCombined = filterC isCombined .|
        (yield header >> mapC edgeToLine) .| unlinesAsciiC .| sinkFile combined
      where
        header = ":START_ID,:END_ID,weight,:TYPE"
    outputBindSites = filterC (not . isCombined) .|
        (yield header >> mapC edgeToLine) .| unlinesAsciiC .| sinkFile binding
      where
        header = ":START_ID,:END_ID,chr,start:int,end:int," <>
            "annotation,affinity,:TYPE"
    isCombined e = case _edge_type e of
        Combined _ -> True
        Binding{..} -> False
{-# INLINE sinkLinkage #-}

createLinkage_ :: Monad m
               => BEDTree Double   -- ^ Promoter activities
               -> BEDTree Double   -- ^ Enhancer activities
               -> M.HashMap GeneName (Double, Double)   -- ^ Gene expression
               -> ConduitT (GeneName, ([BED], [BED]))
                           NetEdge
                           (StateT (S.Set NetNode) m) ()
createLinkage_ act_pro act_enh expr = concatMapMC $ \(geneName, (ps, es)) -> do
    let tfEnhancer = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (getWeight act_enh) es
        edgeEnhancer = flip concatMap tfEnhancer $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "enhancer"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfPromoter = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (getWeight act_pro) ps
        edgePromoter = flip concatMap tfPromoter $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "promoter"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfs = M.toList $ fmap (lp 2 . map (fromJust . (^.score))) $
            M.fromListWith (++) $ tfEnhancer ++ tfPromoter
        (geneExpr, scaledGeneExpr) = M.lookupDefault (0.1, 0) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_expression = Just geneExpr
                           , _node_scaled_expression = Just scaledGeneExpr }
    modify' $ S.insert geneNode
    edgeCombined <- forM tfs $ \(tfName, w) -> do
        let (tfExpr, scaledTfExpr) = M.lookupDefault (0.1, 0) tfName expr
            tfNode = NetNode { _node_name = tfName
                             , _node_expression = Just tfExpr
                             , _node_scaled_expression = Just scaledTfExpr }
        modify' $ S.insert tfNode
        return $ NetEdge { _edge_from = tfName
                         , _edge_to = geneName
                         , _edge_type = Combined (w * sqrt tfExpr) }
    return $ edgePromoter ++ edgeEnhancer ++ edgeCombined
  where
    getBestMotif xs = runIdentity $ runConduit $
        mergeBedWith (maximumBy (comparing (^.score))) xs .| sinkList
    getWeight act bed = case IM.elems (intersecting act bed) of
        [] -> Nothing
        xs -> let p = transform_peak_height $ maximum xs
                  w = sqrt $ transform_site_pvalue (fromJust $ bed^.score) * p
              in if w < edge_weight_cutoff
                then Nothing else Just (getTFName bed, [bed & score .~ Just w])
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE createLinkage_ #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}

-- | Build the network from files containing the information of nodes and edges.
mkNetwork :: FilePath  -- ^ nodes
          -> FilePath  -- ^ edges
          -> IO (Graph 'D NetNode Double)
mkNetwork nodeFl edgeFl = do
    nodeMap <- M.fromList . map ((_node_name &&& id) . nodeFromLine) .
        tail . B.lines <$> B.readFile nodeFl
    runResourceT $ fromLabeledEdges' edgeFl (toEdge nodeMap)
  where
    toEdge nodeMap fl = sourceFileBS fl .| linesUnboundedAsciiC .|
        (dropC 1 >> mapC f) 
      where
        f l = ( ( M.lookupDefault undefined (mk f2) nodeMap
                , M.lookupDefault undefined (mk f1) nodeMap )
              , readDouble f3 )
          where
            [f1,f2,f3,_] = B.split ',' l
{-# INLINE mkNetwork #-}

-- | Read network files as nodes and edges
readNodesAndEdges :: FilePath   -- ^ nodes
                  -> FilePath   -- ^ edges
                  -> IO ([NetNode], [((GeneName, GeneName), Double)])
readNodesAndEdges nodeFl edgeFl = do
    nds <- map nodeFromLine . tail . B.lines <$> B.readFile nodeFl
    es <- map f . tail . B.lines <$> B.readFile edgeFl
    return (nds, es)
  where
    f l = ( ( mk f2, mk f1), readDouble f3 )
        where
        [f1,f2,f3,_] = B.split ',' l