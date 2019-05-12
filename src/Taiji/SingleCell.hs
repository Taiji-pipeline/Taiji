{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeFamilies #-}
module Taiji.SingleCell
    ( prepDataSet
    , computeRanksSC) where

import           Bio.Utils.Misc                    (readDouble, readInt)
import           Bio.Data.Experiment
import Control.Monad.State.Strict
import           Bio.Data.Bed
import Data.BBI.BigBed
import           Conduit
import Data.Conduit.Internal (zipSinks)
import           Control.Lens                      hiding (pre, to)
import qualified Data.Set as S
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import qualified Data.HashMap.Strict                   as M
import           Data.Maybe                        (fromJust, catMaybes)
import           Scientific.Workflow               hiding (_data)
import System.IO.Temp (withTempDirectory)

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils hiding (readPromoters)
import Taiji.Pipeline.SC.ATACSeq.Types

import Taiji.Core.Utils
import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.Types

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
getRanks promoters expr tags bbidx = withTempDirectory "./" "tmp_rank" $ \dir -> do
    let edgeFl = dir ++ "/edge_file"
        nodeFl = dir ++ "/node_file" 
    tfbs <- runConduit $ mapM_ (\x -> queryBB x bbidx) tags .| sinkTFBS
    let proc = return () .| findTargets tfbs promoters .|
            createLinkage_ expr .|
            zipSinks (outputCombinedEdges edgeFl) (return ())
    runResourceT (execStateT (runConduit proc) S.empty) >>=
        outputNodes nodeFl
    pageRank' <$> readNetwork nodeFl edgeFl
{-# INLINE getRanks #-}

cutSite2Bed :: CutSite -> BED3
cutSite2Bed (CutSite chr pos) = asBed chr (pos-50) (pos+50) :: BED3

type BBIndex = M.HashMap B.ByteString BBedFile

-- | Open bigbed files.
openBBs :: [(B.ByteString, Maybe (File '[] 'BigBed))]
        -> IO BBIndex
openBBs xs = fmap (M.fromList . catMaybes) $ forM xs $ \(chr, x) -> case x of
    Nothing -> return Nothing
    Just fl -> do
        bb <- openBBedFile $ fl^.location
        return $ Just (chr, bb)

queryBB :: CutSite -> BBIndex -> ConduitT () TFBS IO ()
queryBB site idx = case M.lookup (bed^.chrom) idx of
    Nothing -> return ()
    Just i -> query (bed^.chrom, bed^.chromStart, bed^.chromEnd) i .| mapC f
  where
    bed = cutSite2Bed site
    f (chr, s, e, rest) = BEDExt (asBed chr s e) info
      where
        info = SiteInfo
            { _tf_name = mk $ head $ B.split '+' f1
            , _site_affinity = toSiteAffinity $ readInt f2
            , _peak_affinity = toPeakAffinity 100 }
        (f1:f2:_) = B.split '\t' rest

sinkTFBS :: Monad m => ConduitT TFBS o m (BEDTree [SiteInfo])
sinkTFBS = mapC (\x -> (x^._bed, [x^._data])) .| sinkList >>=
    return . (fmap . fmap) nub' . bedToTree (++)
  where
    nub' = M.elems . M.fromListWith (\a b -> if g a > g b then a else b) .
        map (\x -> (_tf_name x, x))
      where
        g = getSiteAffinity . _site_affinity