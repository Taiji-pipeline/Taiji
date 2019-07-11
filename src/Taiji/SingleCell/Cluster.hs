-- | Taiji analysis at the cluster level
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeFamilies #-}

module Taiji.SingleCell.Cluster where

import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Conduit
import           Control.Monad.Reader              (asks, ReaderT)
import qualified Data.ByteString.Char8             as B
import qualified Data.Text as T
import qualified Data.HashMap.Strict                   as M
import Bio.Pipeline

import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.SingleCell.Utils
import           Taiji.Prelude
import Taiji.Utils

prepareData :: (a  -- ^ TFBS
               , [(B.ByteString, File t1 'NarrowPeak)]  -- ^ Peaks
               , Maybe (File t2 'Tsv) ) -- ^ Expression
            -> [(a, ((B.ByteString, File t1 'NarrowPeak), File t2 'Tsv))]
prepareData (tfFl, peaks, Just expr) = zip (repeat tfFl) $ zip peaks $ repeat expr
prepareData _ = []

computeRanksCluster :: ( [(B.ByteString, Maybe (File '[] 'BigBed))]  -- ^ TFBS
                       , ((B.ByteString, File t1 'NarrowPeak), File t2 'Tsv) )
                    -> ReaderT TaijiConfig IO
                        (T.Text, [(GeneName, (Double, Double))])
computeRanksCluster (tfFl, ((nm, peakFl), expr)) = do
    promoters <- fromJust <$> asks _taiji_annotation >>= liftIO . readPromoters
    liftIO $ do
        expr' <- (fmap . fmap) (\(a,b) -> (logBase 2 $ a + 1, exp b)) $
            readExpression 1 nm $ expr^.location
        openRegions <- readBed $ peakFl ^.location
        idx <- openBBs tfFl
        ranks <- getRanks (findActivePromoters openRegions promoters)
            expr' openRegions idx
        return (T.pack $ B.unpack nm, ranks)

outputRanksCluster :: [(T.Text, [(GeneName, (Double, Double))])]
                   -> ReaderT TaijiConfig IO ()
outputRanksCluster input = do
    dir <- asks _taiji_output_dir >>= getPath . (<> "/Rank_Cluster/")
    let output1 = dir <> "GeneRanks.tsv"
        output2 = dir <> "GeneRanks_PValues.tsv"
        output3 = dir <> "GeneRanks.html"
    liftIO $ outputRanks output1 output2 output3 input

getRanks :: [Promoter]   -- ^ Active promoters
         -> M.HashMap GeneName (Double, Double)   -- ^ Gene expression
         -> [NarrowPeak]   -- ^ ATAC reads
         -> BBIndex   -- ^ TFBS
         -> IO [(GeneName, (Double, Double))]
getRanks promoters expr tags bbidx = do
    tfbs <- runConduit $ mapM_ (\x -> queryBB x bbidx) tags .| mkTFBSMap
    net <- runResourceT $ runConduit $ mempty .| findTargets tfbs promoters .|
        mkNetwork expr 
    pageRank net
{-# INLINE getRanks #-}

