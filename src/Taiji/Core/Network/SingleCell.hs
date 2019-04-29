{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeFamilies #-}
module Taiji.Core.Network.SingleCell
    ( prepDataSet
    , computeRanksSC) where

import           Bio.Utils.Misc                    (readDouble)
import           Bio.Data.Experiment
import Control.Monad.State.Strict
import           Bio.Data.Bed
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
import System.IO.Temp (withTempFile)
import           Data.Singletons.Prelude.List   (Elem)

import Taiji.Pipeline.SC.ATACSeq.Functions.Utils
import Taiji.Pipeline.SC.ATACSeq.Types

import Taiji.Core.Utils
import Taiji.Core.RegulatoryElement
import Taiji.Core.Network.DeNovo
import Taiji.Core.Ranking
import           Taiji.Types

prepDataSet :: ( f
               , [SCATACSeq S (File t2 'Other)]
               , [RNASeq S (File t3 'Tsv, a, b)] )
            -> IO [(f, File t2 'Other, File t3 'Tsv, [B.ByteString])]
prepDataSet (tfbs, atac, rna) = fmap catMaybes $ forM atac $ \e ->
    case M.lookup (fromJust $ e^.groupName) rnaMap of
        Nothing -> return Nothing
        Just fl -> do
            let x = e^.replicates._2.files
            cells <- withCutSiteIndex (x^.location) $ return . getKeys
            return $ Just (tfbs, x, fl, cells)
  where
    rnaMap = M.fromList $
        map (\x -> (fromJust $ x^.groupName, x^.replicates._2.files._1)) rna

type RankResult = [(GeneName, (Double, Double))]

computeRanksSC :: Elem 'Gzip t1 ~ 'True
               => ( [File t1 'Bed]     -- ^ TFBS
                  , File t2 'Other   -- ^ CutSiteIndex
                  , File t3 'Tsv     -- ^ Expression
                  , [B.ByteString] )   -- ^ Cell Barcode
               -> WorkflowConfig TaijiConfig [RankResult]
computeRanksSC (tfFl, idxFl, rna, cells) = do
    promoters <- fromJust <$> asks _taiji_annotation >>= liftIO . readPromoters
    liftIO $ fmap catMaybes $ forM cells $ \cellBc -> do
        sites <- withCutSiteIndex (idxFl^.location) (lookupIndex cellBc)
        expr <- lookupExpr cellBc $ rna^.location
        case (sites, expr) of
            (Just s, Just e) -> Just <$> getRanks 
                (findActivePromoters (map cutSite2Bed s) promoters) e s tfFl
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

getRanks :: Elem 'Gzip tags ~ 'True
         => [Promoter]   -- ^ Active promoters
         -> M.HashMap GeneName (Double, Double)   -- ^ Gene expression
         -> [CutSite]   -- ^ ATAC reads
         -> [File tags 'Bed]   -- ^ TFBS
         -> IO [(GeneName, (Double, Double))]
getRanks promoters expr tags tfbs_fl = withTempFile "./" "tmp_edge_file" $
    \edgeFl _ -> withTempFile "./" "tmp_node_file" $ \nodeFl _ -> do
        tfbs <- runResourceT $ runConduit $
            mapM_ (streamBedGzip . (^.location)) tfbs_fl .|
            getTFBS (fromCutSites tags)
        let proc = return () .| findTargets tfbs promoters .|
                createLinkage_ expr .|
                zipSinks (outputCombinedEdges edgeFl) (return ())
        runResourceT (execStateT (runConduit proc) S.empty) >>=
            outputNodes nodeFl
        readNetwork nodeFl edgeFl >>= pageRank
{-# INLINE getRanks #-}

fromCutSites :: [CutSite] -> BEDTree PeakAffinity
fromCutSites = bedToTree max . map f
  where
    f x = (cutSite2Bed x, sc)
    sc = toPeakAffinity 100
{-# INLINE fromCutSites #-}

cutSite2Bed :: CutSite -> BED3
cutSite2Bed (CutSite chr pos) = asBed chr (pos-50) (pos+50) :: BED3