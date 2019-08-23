{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Ranking
    ( computeRanks
    , outputRanks
    , pageRank
    , pageRank'
    ) where

import           Control.Monad
import           Control.Monad.ST (runST)
import           Control.Monad.Reader              (asks, ReaderT)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (original)
import           Data.List                         (sort)
import           Data.List.Ordered                 (nubSort)
import qualified Data.HashMap.Strict                   as M
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import qualified Data.Vector as V
import           IGraph
import           IGraph.Algorithms (pagerank, rewireEdges)
import IGraph.Random (withSeed)
import Data.Vector.Algorithms.Search (binarySearch)
import Bio.Utils.Functions (scale)
import Statistics.Sample

import Taiji.Utils.Plot
import Taiji.Utils.Plot.ECharts
import qualified Taiji.Utils.DataFrame as DF
import           Taiji.Core.Network
import           Taiji.Utils
import           Taiji.Prelude

computeRanks :: ( ATACSeq S (File '[] 'Other, File '[] 'Other)
                , Maybe (File '[] 'Tsv) )        -- ^ Expression
             -> ReaderT TaijiConfig IO (T.Text, [(GeneName, (Double, Double))])
computeRanks (atac, exprFl) = do
    maybeNet <- asks _taiji_external_network
    gr <- liftIO $ case maybeNet of
        Nothing -> readNetwork (nodeFl^.location) $ edgeFl^.location
        Just net -> do
            expr <- fromMaybe (return M.empty) $ fmap
                (readExpression 1 (B.pack $ T.unpack grp) . (^.location)) exprFl
            n1 <- getExternalLinks expr net
            n2 <- readAssociations (nodeFl^.location) (edgeFl^.location)
            return $ combineNetwork n1 n2
    result <- liftIO $ pageRank gr
    return (grp, result)
  where
    (nodeFl, edgeFl) = atac^.replicates._2.files
    grp = atac^.groupName._Just

pageRank :: Graph 'D NetNode Double
         -> IO [(GeneName, (Double, Double))]
pageRank gr = do
    getP <- getRankPvalue 10 gr
    return $ flip mapMaybe (zip [0..] (pageRank_ gr)) $ \(i, r) -> if (null $ pre gr i)
        then Nothing
        else Just (_node_name $ nodeLab gr i, (r, getP r))
{-# INLINE pageRank #-}

-- | Do not compute p-values
pageRank' :: Graph 'D NetNode Double -> [(GeneName, Double)]
pageRank' gr = flip mapMaybe (zip [0..] (pageRank_ gr)) $ \(i, r) ->
    if (null $ pre gr i)
        then Nothing
        else Just (_node_name $ nodeLab gr i, r)
{-# INLINE pageRank' #-}

pageRank_ :: Graph 'D NetNode Double
          -> [Double]
pageRank_ gr = pagerank gr 0.85 (Just _node_weight) (Just id)

-- | Compute p-values for PageRank scores, by randomizing the networks.
getRankPvalue :: Int   -- ^ The number of randomization to be performed
              -> Graph 'D NetNode Double
              -> IO (Double -> Double)
getRankPvalue n gr = do
    gr' <- thaw gr
    scores <- fmap (U.fromList . sort . concat) $ replicateM n $
        withSeed 2394 $ \gen -> do
            rewireEdges gr' 1 False False gen
            pageRank_ <$> unsafeFreeze gr' 
    return $ getP scores
  where
    getP vec x = 1 - fromIntegral (bisect vec x) / fromIntegral (U.length vec) 
    bisect v x = runST $ do
        v' <- U.unsafeThaw v
        binarySearch v' x
{-# INLINE getRankPvalue #-}

outputRanks :: FilePath
            -> FilePath
            -> FilePath
            -> [(T.Text, [(GeneName, (Double, Double))])]
            -> IO ()
outputRanks _ _ _ [] = return ()
outputRanks rankFl pValueFl pltFl inputs = do
    DF.writeTable rankFl (T.pack . show) $ fst $ DF.unzip df
    DF.writeTable pValueFl (T.pack . show) $ snd $ DF.unzip df
    savePlots pltFl [] [plt <> toolbox]
  where
    ranks = map (M.fromList . snd) inputs
    genes = nubSort $ concatMap M.keys ranks
    df = DF.transpose $ DF.mkDataFrame (map fst inputs)
        (map (T.pack . B.unpack . original) genes) $ flip map ranks $ \r ->
        flip map genes $ \g -> M.lookupDefault (0,1) g r
    plt = heatmap $ DF.orderDataFrame id $
        DF.mapRows scale $ DF.filterRows (const $ filtCV 0.2) $
        DF.filterRows (const $ V.any (>=1e-5)) $ fst $ DF.unzip df

-- | Determine whether the input pass the CV cutoff
filtCV :: Double -> V.Vector Double -> Bool
filtCV thres xs = sqrt v / m >= thres
  where
    (m, v) = meanVarianceUnb xs