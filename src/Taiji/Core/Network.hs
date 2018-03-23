{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network where

import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.Utils.Functions                             (scale)
import           Control.Monad.Reader                            (asks)
import           Data.Monoid                                     ((<>))
import           Bio.Utils.Misc                                  (readDouble)
import           Bio.Pipeline.Utils                              (asDir,
                                                                  getPath)
import           Conduit
import           Control.Lens                      hiding (pre, to)
import           Control.Monad
import           Data.Binary                       (decodeFile, encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.List                         (transpose)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe                        (fromJust, fromMaybe,
                                                    mapMaybe)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Structure                  (personalizedPagerank)
import           Scientific.Workflow               hiding (_data)

import           Taiji.Core.Config                 ()
import           Taiji.Core.RegulatoryElement      (findTargets)
import           Taiji.Types

-- | Gene and its regulators
type GeneName = CI B.ByteString
type Linkage = (GeneName, [(GeneName, [TFBS])])

computeRanks :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                          , File tag2 'Bed )        -- ^ TFBS
             -> File tag3 'NarrowPeak
             -> File tag4 'NarrowPeak
             -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
             -> Maybe (File '[] 'Tsv)           -- ^ Expression
             -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
computeRanks atac activityPro activityEnh hic expr = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        network <- mkNetwork . createLinkage <$>
            findTargets atac activityPro activityEnh hic

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

createLinkage :: [TFBS] -> [Linkage]
createLinkage tfbs = M.toList $ fmap
    (M.toList . M.fromListWith (++) . map (\x -> (_tf_name $ _ext_data x, [x]))) $
    M.fromListWith (++) $ mapMaybe f tfbs
  where
    f x = let gene = _target_gene $ _ext_data x
          in case gene of
              Nothing -> Nothing
              Just g  -> Just (g, [x])
{-# INLINE createLinkage #-}

mkNetwork :: [Linkage] -> LGraph D NetNode NetEdge
mkNetwork links = fromLabeledEdges $ concatMap toEdge links
  where
    toEdge (gene, tfs) = zipWith f (repeat gene) tfs
    f a b = ( (NetNode a Nothing Nothing Nothing, NetNode (fst b) Nothing Nothing Nothing)
            , NetEdge Nothing Nothing (snd b) )
{-# INLINE mkNetwork #-}

assignWeights :: M.Map (CI B.ByteString) (Double, Double)   -- ^ Gene expression
              -> LGraph D NetNode NetEdge
              -> LGraph D NetNode NetEdge
assignWeights weights gr = mapEdges assignEdgeWeight $
    mapNodes assignNodeWeight gr
  where
    assignNodeWeight _ x =
        let (raw, scaled) = M.findWithDefault (0.1, -10) (nodeName x) weights
        in x{nodeExpression=Just raw, nodeScaledExpression=Just scaled}
    assignEdgeWeight (_, to) x =
        let e = M.findWithDefault (0.1, -10) (nodeName $ nodeLab gr to) weights
        in x{weightExpression=Just $ fst e}
{-# INLINE assignWeights #-}

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
        getSiteWeight sites
{-# INLINE pageRank #-}


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

--------------------------------------------------------------------------------
-- Transformation functions
--------------------------------------------------------------------------------

getSiteWeight :: [TFBS] -> Double
getSiteWeight = maximum . map (\x -> transform_peak_height x * transform_site_pvalue x)
{-# INLINE getSiteWeight #-}

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

-- | Transform peak's p-value into probability score
transform_peak_height :: TFBS -> Double
transform_peak_height tfbs = 1 / (1 + exp (-(x - 5)))
  where
    x = fromJust $ (^._data.peak_signal) tfbs
{-# INLINE transform_peak_height #-}

transform_site_pvalue :: TFBS -> Double
transform_site_pvalue tfbs = 1 / (1 + exp (-(x - 5)))
  where
    x = negate $ logBase 10 $ fromJust $ (^._bed.score) tfbs
{-# INLINE transform_site_pvalue #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

transform_node_weight :: Double -> Double
transform_node_weight = exp
{-# INLINE transform_node_weight #-}
