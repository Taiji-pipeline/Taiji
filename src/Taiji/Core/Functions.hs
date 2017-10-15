{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions
    ( getActivePromoter
    , linkGeneToTFs
    , getTFRanksPrep
    , getTFRanks
    , outputRank
    , buildNet
    ) where

import           Bio.Data.Bed                      (BED (..), BED3 (..),
                                                    BEDLike (..))
import qualified Bio.Data.Bed                      as Bed
import           Bio.Data.Experiment
import           Bio.GO.GREAT                      (AssocRule (..),
                                                    get3DRegulatoryDomains,
                                                    getRegulatoryDomains)
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Instances            ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils                (asDir, getPath)
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble, readInt)
import           Conduit
import           Control.Arrow
import           Control.Lens                      hiding (pre)
import           Control.Monad.Reader              (asks)
import           Data.Binary                       (Binary (..), decodeFile,
                                                    encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Default                      (def)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Function                     (on)
import qualified Data.IntervalMap.Strict           as IM
import           Data.List                         (foldl1', groupBy, sortBy,
                                                    transpose)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import           Data.Maybe                        (fromJust, mapMaybe, fromMaybe)
import           Data.Monoid                       ((<>))
import           Data.Ord                          (comparing)
import qualified Data.Set                          as S
import           Data.Singletons                   (SingI)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Structure                  (pagerank,
                                                    personalizedPagerank)
import           Scientific.Workflow

import           Taiji.Core.Config                 ()
import           Taiji.Pipeline.ATACSeq.Config     (ATACSeqConfig (..))
import           Taiji.Types

type GeneName = CI B.ByteString

instance Binary (CI B.ByteString) where
    put = put . original
    get = fmap mk get

-- | Gene and its regulators
type Linkage = (GeneName, [(GeneName, [BED])])

{-
-- | Call permissive peaks using loose cutoff, aiming for high sensitivity.
callLenientPeak :: SingI tags
                => ATACSeq S (File tags 'Bed)
                -> WorkflowConfig TaijiConfig (ATACSeq S (File '[] 'NarrowPeak))
callLenientPeak input = do
    dir <- asks _atacseq_output_dir >>= getPath . (<> (asDir "/Peaks"))
    let fn output fl = callPeaks output fl Nothing $
            def & cutoff .~ PValue 0.01
                & mode .~ NoModel (-100) 200
    liftIO $ mapFileWithDefName (dir++"/") ".narrowPeak" fn input
    -}

getActivePromoter :: SingI tags
                  => ATACSeq S (File tags 'NarrowPeak)
                  -> WorkflowConfig TaijiConfig (ATACSeq S (File tags 'Bed))
getActivePromoter input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks (asDir . _taiji_output_dir) >>= getPath
    let fun output fl = liftIO $ do
            peaks <- Bed.readBed' $ fl^.location :: IO [BED3]
            tss <- getActiveTSS anno peaks
            Bed.writeBed' output tss
            return $ location .~ output $ emptyFile
    mapFileWithDefName (dir ++ "/") "_active_TSS.bed" fun input

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

linkGeneToTFs :: ATACSeq S ( File tags2 'Bed          -- ^ Active promoters
                           , File tags3 'Bed )        -- ^ TFBS
              -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
linkGeneToTFs atac = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        tfSites <- Bed.readBed' $ tfbs^.location
        activePro <- Bed.readBed' $ activeProFl^.location
        regulators <- runResourceT $ findRegulators activePro Nothing tfSites
        let result = flip map regulators $ \(geneName, tfs) ->
                ( mk geneName
                , map ((head *** id) . unzip) $
                    groupBy ((==) `on` fst) $ sortBy (comparing fst) $
                    map (getTFName &&& id) tfs
                )
            output = dir ++ "/" ++ T.unpack (fromJust $ atac^.groupName) ++
                "_gene_TF_assign.bin"
        encodeFile output (result :: [Linkage])

        printEdgeList (dir ++ "/" ++ T.unpack (fromJust $ atac^.groupName) ++ "_network.tsv")
            result

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

printEdgeList :: FilePath -> [Linkage] -> IO ()
printEdgeList output links = B.writeFile output $ B.unlines $
    flip map links $ \(a,b) -> B.intercalate "\t"
    [original a, B.intercalate "," $ map original $ fst $ unzip b]

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

getTFRanksPrep :: (Maybe (File '[] 'Tsv), [(T.Text, File '[] 'Other)])
               -> IO [( Maybe (M.Map GeneName (Double, Double))
                      , (T.Text, File '[] 'Other) )]
getTFRanksPrep (geneExpr, links) = case geneExpr of
    Nothing -> return $ map (\x -> (Nothing, x)) links
    Just e -> do
        rnaseqData <- readExpression (e^.location) 1
        return $ flip map links $ \(grp, x) ->
            (lookup (B.pack $ T.unpack grp) rnaseqData, (grp, x))

getTFRanks :: ( Maybe (M.Map GeneName (Double, Double))
              , (T.Text, File '[] 'Other) )
           -> IO (T.Text, [(GeneName, Double)])
getTFRanks (expr, (grp, fl)) = do
    links <- decodeFile $ fl^.location
    return (grp, pageRank $ buildNet links expr)

outputRank :: [(T.Text, [(GeneName, Double)])]
           -> WorkflowConfig TaijiConfig FilePath
outputRank results = do
    dir <- asks _taiji_output_dir >>= getPath . asDir
    let output = dir ++ "/GeneRanks.tsv"
    liftIO $ B.writeFile output $ B.unlines $
        header : zipWith toBS (map original genes) (transpose ranks)
    return output
  where
    genes = nubSort $ concatMap (fst . unzip) $ snd $ unzip results
    (groupNames, ranks) = unzip $ flip map results $ \(name, xs) ->
        let geneRanks = M.fromList xs
        in (name, flip map genes $ \g -> M.findWithDefault 0 g geneRanks)
    header = B.pack $ T.unpack $ T.intercalate "\t" $ "Gene" : groupNames
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs

pageRank :: LGraph D (GeneName, Double) Double
         -> [(GeneName, Double)]
pageRank gr = flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
    if i `S.member` tfs
        then Just (fst $ nodeLab gr i, rank)
        else Nothing
  where
    labs = map (nodeLab gr) $ nodes gr
    tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
    ranks = let nodeWeights = map snd labs
                edgeWeights = map (edgeLab gr) $ edges gr
            in personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
{-# INLINE pageRank #-}

buildNet :: [Linkage]
         -> Maybe (M.Map GeneName (Double, Double))  -- ^ absolute value and z-score
         -> LGraph D (GeneName, Double) Double
buildNet links expr = fromLabeledEdges $ concatMap toEdge links
  where
    toEdge (gene, tfs) = zipWith fun (repeat (gene, getNodeWeight gene)) tfs
      where
        fun a (b, beds) =
            let w1 = case expr of
                    Nothing -> 1
                    Just expr' -> sqrt $ fromMaybe 0.01 $ fmap fst $ M.lookup b expr'
                w2 = maximum $ map (fromJust . bedScore) beds
            in ((a, (b, getNodeWeight b)), w1 * w2)
        getNodeWeight x = case expr of
            Nothing -> 1
            Just expr' -> exp $ fromMaybe (-10) $ fmap snd $ M.lookup x expr'
{-# INLINE buildNet #-}

-- | Read RNA expression data
readExpression :: FilePath
               -> Double    -- ^ Threshold to call a gene as non-expressed
               -> IO [( B.ByteString    -- ^ cell type
                     , M.Map GeneName (Double, Double)  -- ^ absolute value and z-score
                     )]
readExpression fl cutoff = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (map readDouble . tail) dat
    return $ zipWith (\a b -> (a, M.fromList $ zip rowNames b)) header $
        transpose $ zipWith zip dataTable $ map computeZscore dataTable
  where
    computeZscore xs
        | all (<cutoff) xs || all (==head xs) xs = replicate (length xs) (-10)
        | otherwise = U.toList $ scale $ U.fromList xs
{-# INLINE readExpression #-}

{-
geneCorrelation :: FilePath -> IO
geneCorrelation fl = do
    c <- B.readFile fl
    let (_:dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (map readDouble . tail) dat
        -}
