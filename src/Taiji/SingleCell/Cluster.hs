-- | Taiji analysis at the cluster level
{-# LANGUAGE DataKinds         #-}
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

prepareData :: Experiment e
            => (a  -- ^ TFBS
               , [e S (B.ByteString, File t1 'NarrowPeak)]  -- ^ Peaks
               , [e S (File t2 'Tsv)] ) -- ^ Expression
            -> [(a, e S ((B.ByteString, File t1 'NarrowPeak), File t2 'Tsv))]
prepareData (tfFl, peaks, expr) = zip (repeat tfFl) $ concatMap split $
    es & mapped.replicates._2.files %~ f 
  where
    f (xs, [y]) = zip xs $ repeat y
    f _ = error "Check your input"
    es = zipExp (concatMap split (mergeExp peaks)) $
        concatMap split $ mergeExp expr

computeRanksCluster :: Experiment e
                    => ( [(B.ByteString, Maybe (File '[] 'BigBed))]  -- ^ TFBS
                       , e S ((B.ByteString, File t1 'NarrowPeak), File t2 'Tsv) )
                    -> ReaderT TaijiConfig IO
                        (e S (T.Text, [(GeneName, (Double, Double))]))
computeRanksCluster (tfFl, input) = do
    promoters <- fromJust <$> asks _taiji_annotation >>= liftIO . readPromoters
    input & replicates.traverse.files %%~ liftIO . ( \((nm, peakFl), expr) -> do
        expr' <- (fmap . fmap) (\(a,b) -> (logBase 2 $ a + 1, exp b)) $
            readExpression 1 nm $ expr^.location
        openRegions <- readBed $ peakFl ^.location
        idx <- openBBs tfFl
        ranks <- getRanks (findActivePromoters openRegions promoters)
            expr' openRegions idx
        return (T.pack $ B.unpack nm, ranks) )

outputRanksCluster :: Experiment e 
                   => [e S (T.Text, [(GeneName, (Double, Double))])]
                   -> ReaderT TaijiConfig IO ()
outputRanksCluster input = forM_ (concatMap split $ mergeExp input) $ \e -> do
    dir <- asks _taiji_output_dir >>=
        getPath . (<> (asDir $ "/Rank_Cluster/" ++ T.unpack (e^.eid)))
    let output1 = dir ++ "/GeneRanks.tsv"
        output2 = dir ++ "/GeneRanks_PValues.tsv"
    liftIO $ outputRanks output1 output2 $ e^.replicates._2.files

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

