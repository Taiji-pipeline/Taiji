{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions
    ( getActivePromoter
    , linkGeneToTFs
    ) where

import           IGraph
import           IGraph.Structure                  (pagerank,
                                                    personalizedPagerank)
import           Bio.Data.Bed            (BED (..), BED3 (..), BEDLike (..))
import qualified Bio.Data.Bed            as Bed
import           Bio.Data.Experiment
import           Bio.GO.GREAT            (AssocRule (..),
                                          get3DRegulatoryDomains,
                                          getRegulatoryDomains)
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances  ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils      (getPath)
import           Bio.Utils.Misc          (readInt)
import           Conduit
import           Control.Arrow
import           Control.Lens hiding (pre)
import           Control.Monad.Reader    (asks)
import           Data.Binary             (encodeFile)
import qualified Data.ByteString.Char8   as B
import           Data.CaseInsensitive    (mk)
import           Data.Function           (on)
import qualified Data.HashMap.Strict     as M
import qualified Data.Set as S
import qualified Data.IntervalMap.Strict as IM
import           Data.List               (foldl1', groupBy, sortBy)
import           Data.List.Ordered       (nubSort)
import           Data.Maybe              (fromJust, mapMaybe)
import           Data.Ord                (comparing)
import           Data.Singletons         (SingI)
import qualified Data.Text               as T
import           Scientific.Workflow
import           System.IO.Temp          (withTempFile)

import           Taiji.Core.Config
import           Taiji.Core.Types        (Linkage, GeneName)

getActivePromoter :: SingI tags
                  => ATACSeq (File tags 'Bed)
                  -> WorkflowConfig TaijiConfig (ATACSeq (File tags 'Bed))
getActivePromoter input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks _taiji_output_dir >>= getPath
    let fun output fl = liftIO $ withTempFile "./" "tmp_macs2_file." $ \tmp _ -> do
            _ <- callPeaks tmp fl Nothing $ cutoff .= QValue 0.1
            peaks <- Bed.readBed' tmp :: IO [BED3]
            tss <- getActiveTSS anno peaks
            Bed.writeBed' output tss
            return $ location .~ output $ emptyFile

    nameWith dir "_active_TSS.bed" fun input

-- | Identify active genes by overlapping their promoters with activity indicators.
getActiveTSS :: BEDLike b
             => FilePath   -- ^ gencode file in GTF format
             -> [b]        -- ^ feature that is used to determine the activity
                           -- of promoters, e.g., H3K4me3 peaks or ATAC-seq peaks
             -> IO [BED]
getActiveTSS input peaks = do
    c <- B.readFile input
    let promoters = map ( \((chr, i, isForward), geneName) ->
            if isForward
                then BED chr (i-5000) (i+1000) (Just geneName)
                     (Just $ fromIntegral i) (Just True)
                else BED chr (i-1000) (i+5000) (Just geneName)
                     (Just $ fromIntegral i) (Just False)
            ) $ map f $ filter g $ map (B.split '\t') $ B.lines c
    return $ nubSort $ map (\x -> BED (chrom x) (truncate $ fromJust $ bedScore x)
        (truncate $ fromJust $ bedScore x) (bedName x) Nothing (bedStrand x)) $
        runIdentity $ yieldMany promoters =$= Bed.intersectBed peaks $$ sinkList
  where
    f xs = let name = B.filter (/='\"') $ last $ B.split ' ' $ head $
                      filter (B.isInfixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")
{-# INLINE getActiveTSS #-}

linkGeneToTFs :: ATACSeq ( File tags2 'Bed          -- ^ Active promoters
                         , File tags3 'Bed )        -- ^ TFBS
              -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
linkGeneToTFs atac = do
    dir <- asks _taiji_output_dir >>= getPath
    liftIO $ do
        tfSites <- Bed.readBed' $ tfbs^.location
        activePro <- Bed.readBed' $ activeProFl^.location
        regulators <- runResourceT $ findRegulators activePro Nothing tfSites
        let result = flip map regulators $ \(geneName, tfs) ->
                ( mk geneName, map ((head *** id) . unzip) $
                    groupBy ((==) `on` fst) $ sortBy (comparing fst) $
                    map (getTFName &&& id) tfs )
            output = dir ++ "/" ++ T.unpack (fromJust $ atac^.groupName) ++
                "_gene_TF_assign.bin"
        encodeFile output (result :: [Linkage])
        return (fromJust $ atac^.groupName, location .~ output $ emptyFile)
  where
    getTFName = mk . head . B.split '+' . fromJust . bedName
    [(activeProFl, tfbs)] = atac^..replicates.folded.files
    {-
    loops = case hic of
        Nothing -> Nothing
        Just x -> if x^.groupName /= atac^.groupName
            then error "Expect the group names to be the same."
            else let [fl] = x^..replicates.folded.files
                 in Just $ read3DContact $ fl^.location
                 -}

findRegulators :: Monad m
                 => [BED]  -- ^ Genes
                 -> Maybe (Source m (BED3, BED3))  -- ^ 3D contacts
                 -> [BED]  -- ^ TFBS
                 -> m [(B.ByteString, [BED])]
findRegulators activePro contacts tfbs = do
    (assign3D, rest) <- case contacts of
        Just contacts' -> do
            let tfbsTree = Bed.bedToTree (++) $ map (\x -> (x, [x])) tfbs
            assignments <- contacts' =$=
                get3DRegulatoryDomains tss 5000 1000 =$=
                mapC ( \(bed, gene) ->
                    (gene, concat $ IM.elems $ Bed.intersecting tfbsTree bed) ) $$
                foldlC (\m (gene, tfs) ->
                    M.insertWith S.union gene (S.fromList tfs) m) M.empty
            let assigned = foldl1' S.union $ M.elems assignments
                unassigned = filter (not . (`S.member` assigned)) tfbs
            return (assignments, unassigned)
        Nothing -> return (M.empty, tfbs)

    assign2D <- fmap (M.fromListWith S.union) $ yieldMany regDomains =$=
        Bed.intersectBedWith S.fromList rest =$= filterC (not . null . snd) =$=
        mapC (first (fromJust . bedName)) $$ sinkList

    return $ M.toList $ fmap S.toList $ M.unionWith S.union assign2D assign3D
  where
    tss = map (\BED{..} ->
        ((_chrom, _chromStart, fromJust _strand), fromJust _name)) activePro
    regDomains = map ( \(b, x) ->
        BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing ) $
        getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000) tss
{-# INLINE findRegulators #-}

pageRank :: Maybe (M.HashMap GeneName (Double, Double))   -- ^ Expression data
         -> LGraph D GeneName ()
         -> [(GeneName, Double)]
pageRank expr gr = flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
    if i `S.member` tfs
        then Just (nodeLab gr i, rank)
        else Nothing
  where
    labs = map (nodeLab gr) $ nodes gr
    tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
    ranks = case expr of
        Just expr' ->
            let lookupExpr x = M.lookupDefault (0.01,-10) x expr'
                nodeWeights = map (exp . snd . lookupExpr) labs
                edgeWeights = map (sqrt . fst . lookupExpr . nodeLab gr . snd) $
                    edges gr
            in personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
        Nothing -> pagerank gr Nothing 0.85

buildNet :: [Linkage] -> IO (LGraph D GeneName ())
buildNet links = return $ fromLabeledEdges $ flip concatMap links $ \(a, b) ->
    zip (zip (repeat a) $ fst $ unzip b) $ repeat ()
{-# INLINE buildNet #-}
