-- | Cell level analysis
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeFamilies #-}
module Taiji.SingleCell.Cell
    ( prepDataSet
    , computeRanksSC) where

import           Bio.Utils.Misc                    (readDouble)
import Control.Monad.State.Strict
import           Bio.Data.Bed
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import qualified Data.HashMap.Strict                   as M
import qualified Data.IntervalMap.Generic.Strict as IM
import Control.Arrow
import           IGraph
import Data.List.Ordered (nubSort)

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils hiding (readPromoters)
import Taiji.Pipeline.SC.ATACSeq.Types

import Taiji.Utils
import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.SingleCell.Utils
import           Taiji.Prelude

prepDataSet :: MonadIO m
            => ( f
               , [SCATACSeq S (File t2 'Other)]
               , [RNASeq S (File t3 'Tsv, a, b)] )
            -> m [(f, File t2 'Other, Maybe (File t3 'Tsv), [B.ByteString])]
prepDataSet (tfbs, atac, rna) = liftIO $ forM atac $ \e -> do
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
               -> ReaderT TaijiConfig IO [RankResult]
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

--------------------------------------------------------------------------------
--
--------------------------------------------------------------------------------

findTargets' :: Monad m
            => BEDTree [GeneName]   -- ^ TF binding sites
            -> [Promoter] -- ^ A list of promoters
            -> ConduitT (BED3, BED3)  -- Chromatin loops
                   (GeneName, [GeneName]) m ()
findTargets' tfbs promoters = getRegulaDomain promoters .|
    mapC (second (nubSort . concatMap f))
  where
    f region = concatMap snd $ IM.toList $ intersecting tfbs $ region^._bed
{-# INLINE findTargets' #-}

mkNetwork' :: Monad m
          => M.HashMap GeneName Double   -- ^ gene expression
          -> ConduitT (GeneName, [GeneName]) o m (Graph 'D NetNode Double)
mkNetwork' expr = fmap fromLabeledEdges $ concatMapC mkEdges .| sinkList
  where
    mkEdges (geneName, tfs) = flip map tfs $ \tf ->
        let w = logBase 2 $ M.lookupDefault 0 tf expr + 1
            tfNode = NetNode { _node_name = tf
                            , _node_weight = w
                            , _node_expression = Nothing }
        in ((geneNode, tfNode), w)
      where
        geneWeight = logBase 2 $ M.lookupDefault 0 geneName expr + 1
        geneNode = NetNode { _node_name = geneName
                           , _node_weight = geneWeight
                           , _node_expression = Nothing }
{-# INLINE mkNetwork' #-}