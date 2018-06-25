{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.RegulatoryElement
    ( Linkage
    , getHiCLoops
    , findActivePromoters
    , createLinkage
    ) where

import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils      (asDir, getPath)
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc          (readInt, readDouble)
import           Bio.Utils.Functions               (scale)
import           Conduit
import           Control.Lens
import           Control.Monad.Reader    (asks)
import qualified Data.ByteString.Char8   as B
import           Data.CaseInsensitive    (mk)
import           Data.Either             (lefts)
import           Data.Function           (on)
import qualified Data.Map.Strict     as M
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.Unboxed               as MU
import qualified Data.Set as S
import qualified Data.IntervalMap.Strict as IM
import System.IO
import Control.Monad.State.Strict
import           Data.List               
import           Data.List.Ordered       (nubSort)
import           Data.Maybe              (fromJust, isNothing, mapMaybe, fromMaybe)
import           Data.Monoid             ((<>))
import           Data.Ord
import           Data.Singletons         (SingI)
import qualified Data.Text               as T
import qualified Data.Vector             as V
import           Scientific.Workflow     hiding (_data)
import           Text.Printf             (printf)

import           Taiji.Core.Config       ()
import           Taiji.Types

type HiCWithSomeFile = HiC N [Either SomeFile (SomeFile, SomeFile)]

getHiCLoops :: [HiCWithSomeFile] -> [HiC S (File '[ChromosomeLoop] 'Bed)]
getHiCLoops inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $
        filter (\x -> ChromosomeLoop `elem` getFileTags x) $ lefts fls

findActivePromoters :: SingI tags
                    => ATACSeq S (File tags 'NarrowPeak)
                    -> WorkflowConfig TaijiConfig (ATACSeq S (File tags 'Bed))
findActivePromoters input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks ((<> "/Promoters") . _taiji_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_active_promoters.bed" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        peaks <- readBed' $ fl^.location :: IO [BED3]
        pro <- readPromoters anno
        writeBed' output $ findActivePromoters_ peaks pro
        return $ location .~ output $ emptyFile
        )

findActivePromoters_ :: BEDLike b
                     => [b]        -- ^ Promoter activity indicator
                     -> [Promoter]
                     -> [Promoter]
findActivePromoters_ bed pro = runIdentity $ runConduit $ yieldMany pro .|
    intersectBed bed .| sinkList
{-# INLINE findActivePromoters_ #-}

-- | Get a list of potential TSS from GTF file
readPromoters :: FilePath -> IO [Promoter]
readPromoters = (fmap . concatMap) fn . readGenes
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x | geneStrand = BEDExt (asBed geneChrom (max 0 $ x - 5000) (x + 1000)) geneName
            | otherwise = BEDExt (asBed geneChrom (max 0 $ x - 1000) (x + 5000)) geneName
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts
{-# INLINE readPromoters #-}

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
    dir <- asks ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir) >>= getPath
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
        withFile netEdges WriteMode $ \h1 -> withFile bindingEdges WriteMode $ \h2 -> do
            B.hPutStrLn h1 ":START_ID,:END_ID,weight,:TYPE"
            B.hPutStrLn h2 $ ":START_ID,:END_ID,chr,start:int,end:int," <>
                "annotation,affinity,:TYPE"
            let conduit = findTargets active_pro tfbs hic .|
                    createLinkage_ activityPro activityEnh expr' .|
                    mapM_C (liftIO . outputEdge h1 h2)
            s <- execStateT (runConduit conduit) S.empty
            let nodeHeader = "geneName:ID,expression,expressionZScore"
            B.writeFile netNodes $ B.unlines $ (nodeHeader:) $
                map nodeToLine $ S.toList s
    return $ atac & replicates.mapped.files .~
        ( emptyFile & location .~ netNodes
        , emptyFile & location .~ netEdges )
  where
    outputEdge h1 h2 e = B.hPutStrLn hdl $ edgeToLine e
      where
        hdl = case _edge_type e of
            Combined _ -> h1
            Binding{..} -> h2
    fl_pro = either (^.location) (^.location) pro
    fl_enh = either (^.location) (^.location) enh
    (active_pro, tfbs) = atac^.replicates._2.files
    grp = atac^.groupName._Just
{-# INLINE createLinkage #-}

createLinkage_ :: BEDTree Double   -- ^ Promoter activities
               -> BEDTree Double   -- ^ Enhancer activities
               -> M.Map GeneName (Double, Double)   -- ^ Gene expression
               -> ConduitT (GeneName, ([BED], [BED]))
                           NetEdge
                           (StateT (S.Set NetNode) IO) ()
createLinkage_ act_pro act_enh expr = concatMapMC $ \(geneName, (ps, es)) -> do
    let tfEnhancer = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (g act_enh) es
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
            mapMaybe (g act_pro) ps
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
        (geneExpr, scaledGeneExpr) = M.findWithDefault (0.1, -10) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_expression = Just geneExpr
                           , _node_scaled_expression = Just scaledGeneExpr }
    modify' $ S.insert geneNode
    edgeCombined <- forM tfs $ \(tfName, w) -> do
        let (tfExpr, scaledTfExpr) = M.findWithDefault (0.1, -10) tfName expr
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
    g act bed = case IM.elems (intersecting act bed) of
        [] -> Nothing
        xs -> Just ( getTFName bed
            , [bed & score.mapped %~ (transform_peak_height (maximum xs) *) . transform_site_pvalue] )
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE createLinkage_ #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}

findTargets :: MonadIO m
            => File tag1 'Bed         -- ^ Active promoters
            -> File tag2 'Bed         -- ^ TFBS
            -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
            -> ConduitT () (GeneName, ([BED], [BED])) m ()
findTargets fl_pro fl_tfbs hic = do
    tfbs <- liftIO $ bedToTree (++) . map (\x -> (x, [x])) <$> readBed' (fl_tfbs^.location)
    promoters <- liftIO $ readBed' $ fl_pro^.location
    let enhancers2D = getRegulatoryDomains 1000000 promoters
    enhancers3D <- case hic of
        Nothing -> return []
        Just fl -> liftIO $ runResourceT $ runConduit $ read3DContact (fl^.location) .|
            getContactingRegions promoters .| sinkList
    let regions = groupRegions $
            zip (repeat Promoter) (mergeDomains promoters) ++
            zip (repeat Enhancer) (mergeDomains $ enhancers2D ++ enhancers3D)
    yieldMany regions .| findTargets_ tfbs
  where
    groupRegions = groupBy ((==) `on` ((^._data) . snd)) .
        sortBy (comparing ((^._data) . snd))

findTargets_ :: Monad m
             => BEDTree [BED]   -- ^ TF binding sites
             -> ConduitT [(DomainType, RegDomain)] (GeneName, ([BED], [BED])) m ()
findTargets_ tfbs = mapC ( \xs ->
    (snd (head xs) ^. _data, divide $ concat $ mapMaybe f xs) )
  where
    divide xs = let (a, b) = partition ((==Promoter) . fst) xs
                in (snd $ unzip a, snd $ unzip b)
    f (ty, region) = case IM.elems (intersecting tfbs region) of
        [] -> Nothing
        xs -> Just $ zip (repeat ty) $ concat xs
{-# INLINE findTargets_ #-}

-- | Merge 2D and 3D domains
mergeDomains :: [RegDomain] -> [RegDomain]
mergeDomains regions = runIdentity $ runConduit $ mergeBedWith f regions .|
    concatC .| sinkList
  where
    f xs = concatMap mergeFn $ groupBy ((==) `on` (^._data)) $
        sortBy (comparing (^._data)) xs
    mergeFn xs =
        let nm = head xs ^._data
            beds = runIdentity $ runConduit $ mergeBed xs .| sinkList
        in map (\x -> x & _data .~ nm) beds
{-# INLINE mergeDomains #-}

-- | Given a gene list , compute the rulatory domain for each gene
getRegulatoryDomains :: Int             -- ^ Extension length. A good default is 1M.
                     -> [Promoter] -- ^ A list of promoters
                     -> [RegDomain] -- ^ Regulatory domains
getRegulatoryDomains ext genes
    | null genes = error "No gene available for domain assignment!"
    | otherwise = concat $ loop $ [Nothing] ++ map Just basal ++ [Nothing]
  where
    loop (a:b:c:rest) = fn a b c : loop (b:c:rest)
    loop _            = []
    fn left (Just bed) right =
        [ bed & chromStart .~ leftPos & chromEnd .~ s
        , bed & chromStart .~ e & chromEnd .~ rightPos ]
      where
        chr = bed^.chrom
        s = bed^.chromStart
        e = bed^.chromEnd
        leftPos
            | isNothing left || chr /= fromJust left ^. chrom = max (s - ext) 0
            | otherwise = min s $ max (s - ext) $ fromJust left ^. chromEnd
        rightPos
            | isNothing right || chr /= fromJust right ^. chrom = e + ext   -- TODO: bound check
            | otherwise = max e $ min (e + ext) $ fromJust right ^. chromStart
    fn _ _ _ = undefined
    basal = V.toList $ fromSorted $ sortBed genes
{-# INLINE getRegulatoryDomains #-}

-- | Read 3D contacts from a file, where each line contains 6 fields separated
-- by Tabs corresponding to 2 interacting loci. Example:
-- chr1 [TAB] 11 [TAB] 100 [TAB] chr2 [TAB] 23 [TAB] 200
read3DContact :: FilePath -> ConduitT i (BED3, BED3) (ResourceT IO) ()
read3DContact input = sourceFileBS input .| linesUnboundedAsciiC .| mapC f .|
    filterC g
  where
    f x = let (chr1:s1:e1:chr2:s2:e2:_) = B.split '\t' x
          in (asBed chr1 (readInt s1) (readInt e1), asBed chr2 (readInt s2) (readInt e2))
    g (x1,x2) | x1^.chrom /= x2^.chrom = True
              | x1^.chromEnd < x2^.chromStart = x2^.chromStart - x1^.chromEnd > 5000
              | x2^.chromEnd < x1^.chromStart = x1^.chromStart - x2^.chromEnd > 5000
              | otherwise = False
{-# INLINE read3DContact #-}

-- | Retrieve regions that contact with the gene's promoter.
getContactingRegions :: Monad m
                     => [Promoter]
                     -> ConduitT (BED3, BED3) RegDomain m ()
getContactingRegions promoters =
    let basal = bedToTree (++) $ flip map promoters $
            \pro -> (pro^._bed, [pro^._data])
     in concatMapC $ \(locA, locB) -> map (BEDExt locB) (intersect basal locA) ++
            map (BEDExt locA) (intersect basal locB)
  where
    intersect t x = nubSort $ concat $ IM.elems $ intersecting t x
{-# INLINE getContactingRegions #-}

-- | Read RNA expression data
readExpression :: Double    -- ^ Threshold to call a gene as non-expressed
               -> B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.Map GeneName (Double, Double)) -- ^ absolute value and z-score
readExpression cutoff ct fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (U.fromList . map ((+pseudoCount) . readDouble) . tail) dat
        dataTable' = MU.fromRows $ zipWith U.zip dataTable $
            map computeZscore dataTable :: MU.Matrix (Double, Double)
        idx = fromMaybe (error $ "Cell type:" ++ B.unpack ct ++ " not found!") $
            lookup ct $ zip header [0..]
    return $ M.fromList $ zip rowNames $ U.toList $ dataTable' `MU.takeColumn` idx
  where
    pseudoCount = 0.1
    computeZscore xs
        | U.length xs == 1 = U.map log xs
        | U.all (<cutoff) xs = U.replicate (U.length xs) (-10)
        | U.length xs == 2 = let fc = log $ U.head xs / U.last xs
                             in U.fromList [fc, negate fc]
        | U.all (== U.head xs) xs = U.replicate (U.length xs) 0
        | otherwise = scale xs
{-# INLINE readExpression #-}

