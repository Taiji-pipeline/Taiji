{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network where

import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (getPath)
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
import           Control.Lens                      hiding (pre, to)
import           Control.Monad
import           Control.Monad.ST (runST)
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Conduit.Cereal               (conduitGet2)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Either                       (fromRight)
import           Data.List                         (transpose, sort)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe                        (fromMaybe, mapMaybe, fromJust)
import           Data.Monoid                       ((<>))
import           Data.Serialize                    (decode, encode, get)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Algorithms (pagerank, rewireEdges)
import           Scientific.Workflow               hiding (_data)
import Data.Vector.Algorithms.Search (binarySearch)

import           Taiji.Core.Config                 ()
import           Taiji.Types

computeRanks :: ( ATACSeq S (File '[] 'Other)
                , Maybe (File '[] 'Tsv) )          -- ^ Expression
             -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
computeRanks (atac, expr) = do
    dir <- asks ((<> "/Network") . _taiji_output_dir) >>= getPath
    liftIO $ do
        let output = dir ++ "/" ++ T.unpack grp ++ "_network.bin"
        network <- case expr of
            Nothing -> mkNetwork input
            Just e -> do
                expr' <- readExpression 1
                    (B.pack $ T.unpack grp ) $ e^.location
                assignWeights expr' <$> mkNetwork input
        getRankPvalue 10 (pageRank network) >>= B.writeFile output . encode
        return (grp, location .~ output $ emptyFile)
  where
    input = atac^.replicates._2.files.location
    grp = atac^.groupName._Just

mkNetwork :: FilePath -> IO (Graph 'D NetNode NetEdge)
mkNetwork input = runResourceT $ fromLabeledEdges' input toEdge
  where
    toEdge fl = sourceFileBS fl .| conduitGet2 get .| concatMapC f
      where
        f (gene, tfs1, tfs2) = zipWith g (repeat gene) $ M.toList $ M.fromListWith (+) $ tfs1 ++ tfs2
        g a b = ( (defaultNode { nodeName = a }, defaultNode { nodeName = fst b })
                , defaultEdge { sites = snd b } )
{-# INLINE mkNetwork #-}

assignWeights :: M.Map (CI B.ByteString) (Double, Double)   -- ^ Gene expression
              -> Graph 'D NetNode NetEdge
              -> Graph 'D NetNode NetEdge
assignWeights weights gr = emap assignEdgeWeight $
    nmap assignNodeWeight gr
  where
    assignNodeWeight (_, x) =
        let (raw, scaled) = M.findWithDefault (0.1, -10) (nodeName x) weights
        in x{nodeExpression=Just raw, nodeScaledExpression=Just scaled}
    assignEdgeWeight ((_, to), x) =
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
        | U.length xs == 1 = U.map log xs
        | U.all (<cutoff) xs = U.replicate (U.length xs) (-10)
        | U.length xs == 2 = let fc = log $ U.head xs / U.last xs
                             in U.fromList [fc, negate fc]
        | U.all (== U.head xs) xs = U.replicate (U.length xs) 0
        | otherwise = scale xs
{-# INLINE readExpression #-}

pageRank :: Graph 'D NetNode NetEdge
         -> Graph 'D NetNode NetEdge
pageRank gr = nmap (\(i, x) -> x{pageRankScore=Just $ ranks U.! i}) gr
  where
    ranks = U.fromList $ pagerank gr 0.85 (Just getNodeW) (Just getEdgeW)
    getNodeW = transform_node_weight . fromMaybe (-10) . nodeScaledExpression
    getEdgeW NetEdge{..} = transform_exp (fromMaybe 0.1 weightExpression) *
        getSiteWeight sites
{-# INLINE pageRank #-}

-- | Compute p-values for PageRank scores, by randomizing the networks.
getRankPvalue :: Int   -- ^ The number of randomization to be performed
              -> Graph 'D NetNode NetEdge
              -> IO (Graph 'D NetNode NetEdge)
getRankPvalue n gr = do
    gr' <- thaw gr
    scores <- fmap (U.fromList . sort . concat) $ replicateM n $ do
        rewireEdges gr' 1 False False
        map (fromJust . pageRankScore . snd) . labNodes . pageRank <$>
            unsafeFreeze gr' 
    return $
        nmap (\(_, x) -> x {pageRankPvalue = getP scores <$> pageRankScore x}) gr
  where
    getP vec x = 1 - fromIntegral (bisect vec x) / fromIntegral (U.length vec) 
    bisect v x = runST $ do
        v' <- U.unsafeThaw v
        binarySearch v' x
{-# INLINE getRankPvalue #-}

outputRanks :: [(T.Text, File '[] 'Other)]
            -> WorkflowConfig TaijiConfig ()
outputRanks inputs = do
    dir <- asks _taiji_output_dir >>= getPath
    let output1 = dir ++ "/GeneRanks.tsv"
        output2 = dir ++ "/GeneRanks_PValues.tsv"

    ranks <- forM inputs $ \(_, fl) -> liftIO $ do
        gr <- fmap (fromRight (error "decode fail") . decode) $
            B.readFile $ fl^.location :: IO (Graph 'D NetNode NetEdge)
        return $! M.fromList $ flip mapMaybe (nodes gr) $ \i ->
            if null (pre gr i)
                then Nothing
                else do
                    let n = nodeLab gr i
                    x <- pageRankScore n
                    p <- pageRankPvalue n
                    return (nodeName n, (x, p))
    let genes = nubSort $ concatMap M.keys ranks
        header = B.pack $ T.unpack $ T.intercalate "\t" $
            "Gene" : fst (unzip inputs)
        ranks' = flip map ranks $ \r -> flip map genes $
            \g -> M.findWithDefault (0,1) g r

    liftIO $ do
        B.writeFile output1 $ B.unlines $ header :
            zipWith toBS (map original genes) (transpose $ (map.map) fst ranks')
        B.writeFile output2 $ B.unlines $ header :
            zipWith toBS (map original genes) (transpose $ (map.map) snd ranks')
  where
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs

--------------------------------------------------------------------------------
-- Transformation functions
--------------------------------------------------------------------------------

getSiteWeight :: Double -> Double
getSiteWeight = id
{-# INLINE getSiteWeight #-}

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

transform_node_weight :: Double -> Double
transform_node_weight = exp
{-# INLINE transform_node_weight #-}
