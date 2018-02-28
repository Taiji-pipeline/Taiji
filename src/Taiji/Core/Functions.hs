{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions
    ( getActivePromoter
    , linkGeneToTFs
    , getTFRanks
    , outputRank
    , linkageToGraphWithDefLabel
    , transform_peak_height
    , buildNet
    , readExpression
    ) where

import           Bio.Data.Bed                      (BED (..), BED3 (..),
                                                    BEDLike (..))
import qualified Bio.Data.Bed                      as Bed
import           Bio.Data.Experiment
import           Bio.GO.GREAT                      (AssocRule (..),
                                                    get3DRegulatoryDomains,
                                                    getRegulatoryDomains)
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
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Function                     (on)
import           Data.Hashable                     (Hashable)
import qualified Data.IntervalMap.Strict           as IM
import           Data.List                         (foldl1', groupBy, sortBy,
                                                    transpose)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe                        (fromJust, mapMaybe)
import           Data.Monoid                       ((<>))
import           Data.Ord                          (comparing)
import qualified Data.Set                          as S
import           Data.Singletons                   (SingI)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Structure                  (personalizedPagerank)
import           Scientific.Workflow
import           Statistics.Correlation.Kendall    (kendall)

import           Taiji.Core.Config                 ()
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

getTFRanks :: ( Maybe (File '[] 'Tsv)
              , (T.Text, File '[] 'Other) )
           -> WorkflowConfig TaijiConfig (T.Text, [(GeneName, Double)])
getTFRanks (expr, (grp, fl)) = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        rnaseqData <- case expr of
            Nothing -> return Nothing
            Just e  -> Just <$> readExpression (e^.location) 1
        links <- decodeFile $ fl^.location
        let gr = buildNet False (B.pack $ T.unpack grp) links rnaseqData
            output = dir ++ "/" ++ T.unpack grp ++ "_network.tsv"
            header = "TF\tTarget\tweight_expression\tweight_peak\tweight_correlation"
        B.writeFile output $ B.unlines $ header : showEdges gr
        return (grp, pageRank gr)
  where
    showEdges gr = map (B.intercalate "\t" . f) $ edges gr
      where
        f (fr, to) = let (w1, w2, w3) = edgeLab gr (fr, to)
                     in [ original $ fst $ nodeLab gr to
                        , original $ fst $ nodeLab gr fr
                        , toShortest $ transform_exp w1
                        , toShortest $ transform_peak_height w2
                        , toShortest $ transform_corr w3 ]

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

pageRank :: LGraph D (GeneName, Double) (Double, Double, Double)
         -> [(GeneName, Double)]
pageRank gr = flip mapMaybe (zip [0..] ranks) $ \(i, rank) ->
    if i `S.member` tfs
        then Just (fst $ nodeLab gr i, rank)
        else Nothing
  where
    labs = map (nodeLab gr) $ nodes gr
    tfs = S.fromList $ filter (not . null . pre gr) $ nodes gr
    ranks = let nodeWeights = map snd labs
                edgeWeights = map (combine . edgeLab gr) $ edges gr
            in personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
    combine (x, y, z) = transform_exp x * transform_peak_height y * transform_corr z
{-# INLINE pageRank #-}

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

transform_peak_height :: Double -> Double
transform_peak_height x = 1 / (1 + exp (-(x - 5)))
{-# INLINE transform_peak_height #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

linkageToGraphWithDefLabel ::(Hashable v, Read v, Eq v, Show v, Show e)
                           => v -> e -> [Linkage] -> LGraph D (GeneName, v) e
linkageToGraphWithDefLabel v e links = fromLabeledEdges $ concatMap toEdge links
  where
    toEdge (gene, tfs) = zipWith ( \a b -> (((a,v), (fst b,v)), e) ) (repeat gene) tfs
{-# INLINE linkageToGraphWithDefLabel #-}

buildNet :: Bool  -- ^ whether to calculate correlation
         -> B.ByteString  -- ^ cell type
         -> [Linkage]
         -> Maybe ( M.Map GeneName Int
                  , M.Map B.ByteString Int
                  , MU.Matrix (Double, Double)
                  )
         -> LGraph D (GeneName, Double) (Double, Double, Double)
buildNet _ _ links Nothing = linkageToGraphWithDefLabel 1 (1, 1, 1) links
buildNet useCor celltype links (Just (rows, cols, table)) = fromLabeledEdges $
    concatMap toEdge links
  where
    toEdge (gene, tfs) = map fun tfs
      where
        idx_cell = M.findWithDefault undefined celltype cols
        expr = table `MU.takeColumn` idx_cell
        idx_gene = M.lookup gene rows
        expr_gene = fmap (table `MU.takeRow`) idx_gene
        fun (tf, beds) =
            let idx_tf = M.lookup tf rows
                weight1 = case idx_tf of
                    Nothing -> 0.1
                    Just j  -> fst $ expr U.! j
                weight2 = maximum $ map (fromJust . bedScore) beds
                weight3 = if useCor
                    then case expr_gene of
                            Nothing -> 0.4
                            Just v1 -> case fmap (table `MU.takeRow`) idx_tf of
                                Nothing -> 0.4
                                Just v2 -> let cor = kendall $ U.zip v1 v2
                                           in if isNaN cor then 0 else cor
                    else 1
                node_weight_gene = getNodeWeight idx_gene
                node_weight_tf = getNodeWeight idx_tf
            in ( ((gene, node_weight_gene), (tf, node_weight_tf))
               , (weight1, weight2, weight3) )
        getNodeWeight x = case x of
            Nothing -> exp (-10)
            Just i  -> exp $ snd $ expr U.! i
{-# INLINE buildNet #-}

-- | Read RNA expression data
readExpression :: FilePath
               -> Double    -- ^ Threshold to call a gene as non-expressed
               -> IO ( M.Map GeneName Int  -- ^ row
                     , M.Map B.ByteString Int  -- ^ column
                     , MU.Matrix (Double, Double)  -- ^ absolute value and z-score
                     )
readExpression fl cutoff = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (U.fromList . map ((+pseudoCount) . readDouble) . tail) dat
    return ( M.fromList $ zip rowNames [0..]
           , M.fromList $ zip header [0..]
           , MU.fromRows $ zipWith U.zip dataTable $ map computeZscore dataTable
           )
  where
    pseudoCount = 0.1
    computeZscore xs
        | U.all (<cutoff) xs || U.all (== U.head xs) xs = U.replicate (U.length xs) (-10)
        | otherwise = scale xs
{-# INLINE readExpression #-}
