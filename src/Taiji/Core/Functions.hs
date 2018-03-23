{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions
    ( getActivePromoter
    , outputRank
    , computeRanks
    , transform_peak_height
    , readExpression
    , getHiCLoops
    ) where

import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.Pipeline.Instances                          ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils                              (asDir,
                                                                  getPath)
import           Bio.Utils.Functions                             (scale)
import           Bio.Utils.Misc                                  (readDouble)
import           Conduit
import           Control.Lens                                    hiding (pre)
import           Control.Monad
import           Control.Monad.Reader                            (asks)
import           Control.Monad.State.Lazy                        (execStateT,
                                                                  modify)
import           Data.Binary                                     (decodeFile,
                                                                  encodeFile)
import qualified Data.ByteString.Char8                           as B
import           Data.CaseInsensitive                            (CI, mk,
                                                                  original)
import           Data.Double.Conversion.ByteString               (toShortest)
import           Data.Either                                     (lefts)
import           Data.List                                       (transpose)
import           Data.List.Ordered                               (nubSort)
import qualified Data.Map.Strict                                 as M
import qualified Data.Matrix.Unboxed                             as MU
import           Data.Maybe                                      (fromJust,
                                                                  fromMaybe,
                                                                  mapMaybe)
import           Data.Monoid                                     ((<>))
import           Data.Singletons                                 (SingI)
import qualified Data.Text                                       as T
import qualified Data.Vector.Unboxed                             as U
import           IGraph
import           IGraph.Structure                                (personalizedPagerank)
import           Scientific.Workflow hiding (_data)

import           Taiji.Core.Config                               ()
import           Taiji.Core.Functions.Internal.RegulatoryElement
import           Taiji.Types

-- | Gene and its regulators
type GeneName = CI B.ByteString
type Linkage = (GeneName, [(GeneName, [TFBS])])

type HiCWithSomeFile = HiC N [Either SomeFile (SomeFile, SomeFile)]

getHiCLoops :: [HiCWithSomeFile] -> [HiC S (File '[ChromosomeLoop] 'Bed)]
getHiCLoops inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $
        filter (\x -> ChromosomeLoop `elem` getFileTags x) $ lefts fls

getActivePromoter :: SingI tags
                  => ATACSeq S (File tags 'NarrowPeak)
                  -> WorkflowConfig TaijiConfig (ATACSeq S (File tags 'Bed))
getActivePromoter input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks (asDir . _taiji_output_dir) >>= getPath
    let fun output fl = liftIO $ do
            peaks <- readBed' $ fl^.location :: IO [BED3]
            pro <- readPromoters anno
            writeBed' output $ findActivePromoters peaks pro
            return $ location .~ output $ emptyFile
    mapFileWithDefName (dir ++ "/") "_active_promoters.bed" fun input

findTargets :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                         , File tag2 'Bed )        -- ^ TFBS
            -> Maybe (HiC S (File '[ChromosomeLoop] 'Bed))  -- ^ HiC loops
            -> IO [TFBS]
findTargets atac hic = do
    promoters <- readBed' $ fl_pro^.location
    let enhancers = getRegulatoryDomains 1000000 promoters
    distalEhancers <- case hic of
        Nothing -> return []
        Just x -> if x^.groupName /= atac^.groupName
            then error "Expect the group names to be the same."
            else let [fl] = x^..replicates.folded.files
                 in runResourceT $ runConduit $ read3DContact (fl^.location) .|
                        getContactingRegions promoters .| sinkList

    let action = runConduit $ readBed (fl_tfbs^.location) .| mapC toSite .|
            assignTFBS promoters .| concatMapMC (fn Pro) .|
            assignTFBS distalEhancers .| concatMapMC (fn NearbyEnh) .|
            assignTFBS enhancers .| concatMapMC (fn DistalEnh) .| sinkNull
    execStateT action []
  where
    fn _ (Left x) = return $ Just x
    fn t (Right x) = modify ((_data.target_through .~ Just t $ x) :) >> return Nothing
    toSite x = BEDExt (x :: BED) $ SiteInfo (mk $ fromJust $ x^.name)
        Nothing Nothing Nothing
    [(fl_pro, fl_tfbs)] = atac^..replicates.folded.files
{-# INLINE findTargets #-}

computeRanks :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                          , File tag2 'Bed )        -- ^ TFBS
             -> Maybe (HiC S (File '[ChromosomeLoop] 'Bed))  -- ^ HiC loops
             -> Maybe (File '[] 'Tsv)           -- ^ Expression
             -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
computeRanks atac hic expr = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        network <- mkNetwork . createLinkage <$> findTargets atac hic

        gr <- fmap pageRank $ case expr of
            Nothing -> return network
            Just e -> do
                expr' <- readExpression 1
                    (B.pack $ T.unpack $ fromJust $ atac^.groupName) $ e^.location
                return $ assignWeights expr' network

        let output = dir ++ "/" ++ T.unpack (fromJust $ atac^.groupName) ++
                "_network.bin"
        encodeFile output gr
        return (fromJust $ atac^.groupName, location .~ output $ emptyFile)

outputRank :: [(T.Text, File '[] 'Other)]
           -> WorkflowConfig TaijiConfig FilePath
outputRank inputs = do
    dir <- asks _taiji_output_dir >>= getPath . asDir
    let output = dir ++ "/GeneRanks.tsv"

    ranks <- forM inputs $ \(_, fl) -> liftIO $ do
        gr <- decodeFile $ fl^.location :: IO (LGraph D NetNode NetEdge)
        return $! M.fromList $ flip mapMaybe (nodes gr) $ \i ->
            if null (pre gr i)
                then Nothing
                else do
                    let n = nodeLab gr i
                    x <- pageRankScore n
                    return (nodeName n, x)
    let genes = nubSort $ concatMap M.keys ranks
        header = B.pack $ T.unpack $ T.intercalate "\t" $
            "Gene" : fst (unzip inputs)
        ranks' = flip map ranks $ \r -> flip map genes $
            \g -> M.findWithDefault 0 g r

    liftIO $ B.writeFile output $ B.unlines $
        header : zipWith toBS (map original genes) (transpose ranks')
    return output
  where
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs

pageRank :: LGraph D NetNode NetEdge
         -> LGraph D NetNode NetEdge
pageRank gr = mapNodes (\i x -> x{pageRankScore=Just $ ranks U.! i}) gr
  where
    labs = map (nodeLab gr) $ nodes gr
    nodeWeights = map (transform_node_weight . fromMaybe (-10) .
        nodeScaledExpression) labs
    edgeWeights = map (combine . edgeLab gr) $ edges gr
    ranks = U.fromList $ personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
    combine NetEdge{..} = transform_exp (fromMaybe 1 weightExpression) *
        transform_peak_height sites
{-# INLINE pageRank #-}

mkNetwork :: [Linkage] -> LGraph D NetNode NetEdge
mkNetwork links = fromLabeledEdges $ concatMap toEdge links
  where
    toEdge (gene, tfs) = zipWith f (repeat gene) tfs
    f a b = ( (NetNode a Nothing Nothing Nothing, NetNode (fst b) Nothing Nothing Nothing)
            , NetEdge Nothing Nothing (snd b) )

assignWeights :: M.Map (CI B.ByteString) (Double, Double)   -- ^ Gene expression
              -> LGraph D NetNode NetEdge
              -> LGraph D NetNode NetEdge
assignWeights weights gr = mapEdges assignEdgeWeight $
    mapNodes assignNodeWeight gr
  where
    assignNodeWeight _ x =
        let (raw, scaled) = M.findWithDefault (0.1, -10) (nodeName x) weights
        in x{nodeExpression=Just raw, nodeScaledExpression=Just scaled}
    assignEdgeWeight (fr, to) x =
        let e = M.findWithDefault (0.1, -10) (nodeName $ nodeLab gr to) weights
        in x{weightExpression=Just $ fst e}

-- | Read RNA expression data
readExpression :: Double    -- ^ Threshold to call a gene as non-expressed
               -> B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.Map (CI B.ByteString) (Double, Double)) -- ^ absolute value and z-score
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
        | U.all (<cutoff) xs || U.all (== U.head xs) xs = U.replicate (U.length xs) (-10)
        | otherwise = scale xs
{-# INLINE readExpression #-}



--------------------------------------------------------------------------------
-- Transformation functions
--------------------------------------------------------------------------------

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

transform_peak_height :: [TFBS] -> Double
transform_peak_height beds = 1 / (1 + exp (-(x - 5)))
  where
    x = maximum $ map (fromJust . (^.score) . _ext_bed) beds
{-# INLINE transform_peak_height #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

transform_node_weight :: Double -> Double
transform_node_weight = exp
{-# INLINE transform_node_weight #-}
