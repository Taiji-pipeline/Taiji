{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network where

import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (getPath)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
import Control.Arrow ((&&&))
import           Control.Lens                      hiding (pre, to)
import           Control.Monad
import           Control.Monad.ST (runST)
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.List                         (transpose, sort)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import           Data.Maybe                        (fromMaybe, mapMaybe)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Algorithms (pagerank, rewireEdges)
import           Scientific.Workflow               hiding (_data)
import Data.Vector.Algorithms.Search (binarySearch)

import           Taiji.Core.Config                 ()
import           Taiji.Types

computeRanks :: ATACSeq S (File '[] 'Other, File '[] 'Other)
             -> IO (T.Text, [(GeneName, (Double, Double))])
computeRanks atac = do
    gr <- mkNetwork (nodeFl^.location) $ edgeFl^.location
    result <- pageRank gr
    return (grp, result)
  where
    (nodeFl, edgeFl) = atac^.replicates._2.files
    grp = atac^.groupName._Just

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
        f l = ( ( M.findWithDefault undefined (mk f2) nodeMap
                , M.findWithDefault undefined (mk f1) nodeMap )
              , readDouble f3 )
          where
            [f1,f2,f3,_] = B.split ',' l
{-# INLINE mkNetwork #-}

pageRank :: Graph 'D NetNode Double
         -> IO [(GeneName, (Double, Double))]
pageRank gr = do
    getP <- getRankPvalue 10 gr
    return $ flip mapMaybe (zip [0..] (pageRank_ gr)) $ \(i, r) -> if (null $ pre gr i)
        then Nothing
        else Just (_node_name $ nodeLab gr i, (r, getP r))
{-# INLINE pageRank #-}

pageRank_ :: Graph 'D NetNode Double
          -> [Double]
pageRank_ gr = pagerank gr 0.85 (Just getNodeW) (Just id)
  where
    getNodeW = transform_node_weight . fromMaybe (-10) . _node_scaled_expression

-- | Compute p-values for PageRank scores, by randomizing the networks.
getRankPvalue :: Int   -- ^ The number of randomization to be performed
              -> Graph 'D NetNode Double
              -> IO (Double -> Double)
getRankPvalue n gr = do
    gr' <- thaw gr
    scores <- fmap (U.fromList . sort . concat) $ replicateM n $ do
        rewireEdges gr' 1 False False
        pageRank_ <$> unsafeFreeze gr' 
    return $ getP scores
  where
    getP vec x = 1 - fromIntegral (bisect vec x) / fromIntegral (U.length vec) 
    bisect v x = runST $ do
        v' <- U.unsafeThaw v
        binarySearch v' x
{-# INLINE getRankPvalue #-}

outputRanks :: [(T.Text, [(GeneName, (Double, Double))])]
            -> WorkflowConfig TaijiConfig ()
outputRanks inputs = do
    dir <- asks _taiji_output_dir >>= getPath
    let output1 = dir ++ "/GeneRanks.tsv"
        output2 = dir ++ "/GeneRanks_PValues.tsv"
        ranks = map (M.fromList . snd) inputs
        genes = nubSort $ concatMap M.keys ranks
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
