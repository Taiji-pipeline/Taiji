-- | Taiji analysis at the cluster level
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeFamilies #-}

module Taiji.SingleCell.Cluster where

import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Conduit
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import qualified Data.Text as T
import qualified Data.HashMap.Strict                   as M
import Bio.Pipeline

import Taiji.Core.Utils
import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.SingleCell.Utils
import           Taiji.Prelude

computeRanksCluster :: Experiment e
                    => ( [(B.ByteString, Maybe (File '[] 'BigBed))]  -- ^ TFBS
                       , e S [(B.ByteString, File '[] 'NarrowPeak)] )
                    -> WorkflowConfig TaijiConfig ()
computeRanksCluster ([], _) = return ()
computeRanksCluster (tfFl, exps) = do
    promoters <- fromJust <$> asks _taiji_annotation >>= liftIO . readPromoters
    dir <- asks _taiji_output_dir >>=
        getPath . (<> (asDir $ "/Rank_Cluster/" ++ T.unpack (exps^.eid)))
    let output1 = dir ++ "/GeneRanks.tsv" 
        output2 = dir ++ "/GeneRanks_PValues.tsv"
    liftIO $ do
        res <- forM (exps^.replicates._2.files) $ \(nm, peakFl) -> do
            openRegions <- readBed $ peakFl ^.location
            idx <- openBBs tfFl
            ranks <- getRanks (findActivePromoters openRegions promoters)
                M.empty openRegions idx
            return (T.pack $ B.unpack nm, ranks)
        outputRanks output1 output2 res

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

