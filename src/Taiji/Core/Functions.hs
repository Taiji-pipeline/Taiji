{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions
    ( getActivePromoter
    , outputRank
    , computeRanks
    , transform_peak_height
    , readExpression
    , getHiCLoops
    ) where

import           Bio.Data.Bed                      (BED (..), BED3 (..),
                                                    BEDLike (..))
import qualified Bio.Data.Bed                      as Bed
import           Bio.Data.Experiment
import           Bio.GO.GREAT                      (AssocRule (..),
                                                    get3DRegulatoryDomains,
                                                    getRegulatoryDomains,
                                                    read3DContact)
import           Bio.Pipeline.Instances            ()
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils                (asDir, getPath)
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble, readInt)
import           Conduit
import           Control.Arrow
import           Control.Lens                      hiding (pre)
import           Control.Monad
import           Control.Monad.Reader              (asks)
import           Data.Binary                       (Binary (..), decodeFile,
                                                    encodeFile)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Either                       (lefts)
import           Data.Function                     (on)
import           Data.Hashable                     (Hashable)
import qualified Data.IntervalMap.Strict           as IM
import           Data.List                         (foldl1', groupBy, sortBy,
                                                    transpose)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe                        (fromJust, fromMaybe,
                                                    mapMaybe)
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

-- | Gene and its regulators
type GeneName = CI B.ByteString
type Linkage = (GeneName, [(GeneName, [BED])])

type HiCWithSomeFile = HiC N [Either SomeFile (SomeFile, SomeFile)]

getHiCLoops :: [HiCWithSomeFile] -> [HiC S (File '[ChromosomeLoop] 'Bed)]
getHiCLoops inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $
        filter (\x -> ChromosomeLoop `elem` getFileTags x) $ lefts fls

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

computeRanks :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                           , File tag2 'Bed )        -- ^ TFBS
             -> Maybe (File '[] 'Tsv)           -- ^ Expression
             -> Maybe (HiC S (File '[ChromosomeLoop] 'Bed))  -- ^ HiC loops
             -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
computeRanks atac expr hic = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        tfSites <- Bed.readBed' $ tfbs^.location
        activePro <- Bed.readBed' $ activeProFl^.location
        network <- fmap mkNetwork $ runResourceT $
            findRegulators activePro loops tfSites

        gr <- fmap pageRank $ case expr of
            Nothing -> return network
            Just e -> do
                expr' <- readExpression 1 (B.pack $ T.unpack $ fromJust name) $ e^.location
                return $ assignWeights expr' network

        let output = dir ++ "/" ++ T.unpack (fromJust $ atac^.groupName) ++
                "_network.bin"
        encodeFile output gr
        return (fromJust $ atac^.groupName, location .~ output $ emptyFile)
  where
    name = atac^.groupName
    [(activeProFl, tfbs)] = atac^..replicates.folded.files
    loops = case hic of
        Nothing -> Nothing
        Just x -> if x^.groupName /= name
            then error "Expect the group names to be the same."
            else let [fl] = x^..replicates.folded.files
                 in Just $ read3DContact $ fl^.location

findRegulators :: Monad m
               => [BED]  -- ^ Genes
               -> Maybe (Source m (BED3, BED3))  -- ^ 3D contacts
               -> [BED]  -- ^ TFBS
               -> m [Linkage]
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

    let links = M.toList $ fmap S.toList $ M.unionWith S.union assign2D assign3D
    return $ flip map links $ \(geneName, tfs) ->
        ( mk geneName
        , map ((head *** id) . unzip) $ groupBy ((==) `on` fst) $
            sortBy (comparing fst) $ map (getTFName &&& id) tfs
        )
  where
    tss = map (\BED{..} ->
        ((_chrom, _chromStart, fromJust _strand), fromJust _name)) activePro
    regDomains = map ( \(b, x) ->
        BED (chrom b) (chromStart b) (chromEnd b) (Just x) Nothing Nothing ) $
        getRegulatoryDomains (BasalPlusExtension 5000 1000 1000000) tss
    getTFName = mk . head . B.split '+' . fromJust . bedName
{-# INLINE findRegulators #-}

{-
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
                        -}

outputRank :: [(T.Text, File '[] 'Other)]
           -> WorkflowConfig TaijiConfig FilePath
outputRank inputs = do
    dir <- asks _taiji_output_dir >>= getPath . asDir
    let output = dir ++ "/GeneRanks.tsv"

    ranks <- forM inputs $ \(ct, fl) -> liftIO $ do
        gr <- decodeFile $ fl^.location :: IO (LGraph D NetNode NetEdge)
        return $! M.fromList $ flip mapMaybe (nodes gr) $ \i ->
            if null (pre gr i)
                then Nothing
                else do
                    let n = nodeLab gr i
                    x <- pageRankScore n
                    return (nodeName n, x)
    let genes = nubSort $ concatMap M.keys ranks
        header = B.pack $ T.unpack $ T.intercalate "\t" $
            "Gene" : fst (unzip inputs)
        ranks' = flip map ranks $ \r -> flip map genes $
            \g -> M.findWithDefault 0 g r

    liftIO $ B.writeFile output $ B.unlines $
        header : zipWith toBS (map original genes) (transpose ranks')
    return output
  where
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs

pageRank :: LGraph D NetNode NetEdge
         -> LGraph D NetNode NetEdge
pageRank gr = mapNodes (\i x -> x{pageRankScore=Just $ ranks U.! i}) gr
  where
    labs = map (nodeLab gr) $ nodes gr
    nodeWeights = map (transform_node_weight . fromMaybe (-10) .
        nodeScaledExpression) labs
    edgeWeights = map (combine . edgeLab gr) $ edges gr
    ranks = U.fromList $ personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
    combine NetEdge{..} = transform_exp (fromMaybe 1 weightExpression) *
        transform_peak_height sites
{-# INLINE pageRank #-}

mkNetwork :: [Linkage] -> LGraph D NetNode NetEdge
mkNetwork links = fromLabeledEdges $ concatMap toEdge links
  where
    toEdge (gene, tfs) = zipWith f (repeat gene) tfs
    f a b = ( (NetNode a Nothing Nothing Nothing, NetNode (fst b) Nothing Nothing Nothing)
            , NetEdge Nothing Nothing (snd b) )

assignWeights :: M.Map (CI B.ByteString) (Double, Double)
              -> LGraph D NetNode NetEdge
              -> LGraph D NetNode NetEdge
assignWeights weights gr = mapEdges assignEdgeWeight $
    mapNodes assignNodeWeight gr
  where
    assignNodeWeight _ x =
        let (raw, scaled) = M.findWithDefault (0.1, -10) (nodeName x) weights
        in x{nodeExpression=Just raw, nodeScaledExpression=Just scaled}
    assignEdgeWeight (fr, to) x =
        let e = M.findWithDefault (0.1, -10) (nodeName $ nodeLab gr to) weights
        in x{weightExpression=Just $ fst e}


-- | Read RNA expression data
readExpression :: Double    -- ^ Threshold to call a gene as non-expressed
               -> B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.Map (CI B.ByteString) (Double, Double)) -- ^ absolute value and z-score
readExpression cutoff ct fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (U.fromList . map ((+pseudoCount) . readDouble) . tail) dat
        dataTable' = MU.fromRows $ zipWith U.zip dataTable $
            map computeZscore dataTable :: MU.Matrix (Double, Double)
        idx = fromMaybe (error $ "Cell type:" ++ B.unpack ct ++ " not found!") $
            lookup ct $ zip header [0..]
    return $ M.fromList $ zip rowNames $ U.toList $ dataTable' `MU.takeColumn` idx
  where
    pseudoCount = 0.1
    computeZscore xs
        | U.all (<cutoff) xs || U.all (== U.head xs) xs = U.replicate (U.length xs) (-10)
        | otherwise = scale xs
{-# INLINE readExpression #-}



--------------------------------------------------------------------------------
-- Transformation functions
--------------------------------------------------------------------------------

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

transform_peak_height :: [BED] -> Double
transform_peak_height beds = 1 / (1 + exp (-(x - 5)))
  where
    x = maximum $ map (fromJust . bedScore) beds
{-# INLINE transform_peak_height #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

transform_node_weight :: Double -> Double
transform_node_weight = exp
{-# INLINE transform_node_weight #-}
