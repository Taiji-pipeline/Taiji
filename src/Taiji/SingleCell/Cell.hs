-- | Cell level analysis
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeFamilies #-}
module Taiji.SingleCell.Cell
    ( prepDataSet
    , computeRanksSC) where

import           Bio.Utils.Misc                    (readDouble)
import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import qualified Data.HashMap.Strict                   as M
import           Scientific.Workflow               hiding (_data)

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils hiding (readPromoters)
import Taiji.Pipeline.SC.ATACSeq.Types

import Taiji.Core.Utils
import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.SingleCell.Utils
import           Taiji.Prelude

prepDataSet :: ( f
               , [SCATACSeq S (File t2 'Other)]
               , [RNASeq S (File t3 'Tsv, a, b)] )
            -> IO [(f, File t2 'Other, Maybe (File t3 'Tsv), [B.ByteString])]
prepDataSet (tfbs, atac, rna) = forM atac $ \e -> do
    let x = e^.replicates._2.files
    cells <- withCutSiteIndex (x^.location) $ return . getKeys
    return (tfbs, x, M.lookup (fromJust $ e^.groupName) rnaMap, cells)
  where
    rnaMap = M.fromList $
        map (\x -> (fromJust $ x^.groupName, x^.replicates._2.files._1)) rna

type RankResult = [(GeneName, Double)]

computeRanksSC :: ( [(B.ByteString, Maybe (File '[] 'BigBed))]  -- ^ TFBS
                  , File t2 'Other   -- ^ CutSiteIndex
                  , Maybe (File t3 'Tsv)     -- ^ Expression
                  , [B.ByteString] )   -- ^ Cell Barcode
               -> WorkflowConfig TaijiConfig [RankResult]
computeRanksSC (tfFl, idxFl, rna, cells) = do
    promoters <- fromJust <$> asks _taiji_annotation >>= liftIO . readPromoters
    idx <- liftIO $ openBBs tfFl
    liftIO $ fmap catMaybes $ forM cells $ \cellBc -> do
        print cellBc
        sites <- withCutSiteIndex (idxFl^.location) (lookupIndex cellBc)
        print $ fmap length sites
        expr <- maybe (return $ Just M.empty) (lookupExpr cellBc . (^.location)) rna
        case (sites, expr) of
            (Just s, Just e) -> Just <$> getRanks 
                (findActivePromoters (map cutSite2Bed s) promoters) e s idx
            _ -> return Nothing

lookupExpr :: B.ByteString
           -> FilePath
           -> IO (Maybe (M.HashMap GeneName (Double, Double)))
lookupExpr ct input = do
    (header:content) <- B.lines <$> B.readFile input
    let geneNames = map mk $ B.split '\t' header
    return $ do
        r <- lookup ct $ map ((\(x:xs) -> (x, xs)) . B.split '\t') content
        return $ M.fromList $ zip geneNames $
            map (\x -> let x' = logBase 2 $ readDouble x + 1 in (x',x')) r

getRanks :: [Promoter]   -- ^ Active promoters
         -> M.HashMap GeneName (Double, Double)   -- ^ Gene expression
         -> [CutSite]   -- ^ ATAC reads
         -> BBIndex   -- ^ TFBS
         -> IO [(GeneName, Double)]
getRanks promoters expr tags bbidx = do
    tfbs <- runConduit $ mapM_ (\x -> queryBB (cutSite2Bed x) bbidx) tags .|
        mkTFBSMap
    fmap pageRank' $ runResourceT $ runConduit $ mempty .| findTargets tfbs promoters .|
        mkNetwork expr 
{-# INLINE getRanks #-}

cutSite2Bed :: CutSite -> BED3
cutSite2Bed (CutSite chr pos) = asBed chr (pos-50) (pos+50) :: BED3