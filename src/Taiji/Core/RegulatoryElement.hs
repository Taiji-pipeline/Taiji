{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.RegulatoryElement
    ( Linkage
    , getHiCLoops
    , findActivePromoters
    , findTargets
    , createLinkage
    ) where

import           Bio.Data.Bed
import           Bio.Data.Experiment
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils      (asDir, getPath)
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc          (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.Reader    (asks)
import qualified Data.ByteString.Char8   as B
import           Data.CaseInsensitive    (mk)
import           Data.Conduit.Cereal     (conduitGet2, conduitPut)
import           Data.Either             (lefts)
import           Data.Function           (on)
import qualified Data.HashMap.Strict     as M
import qualified Data.IntervalMap.Strict as IM
import           Data.List               (foldl', groupBy, partition, sortBy)
import           Data.List.Ordered       (nubSort)
import           Data.Maybe              (fromJust, isNothing, mapMaybe)
import           Data.Monoid             ((<>))
import           Data.Ord
import           Data.Serialize          (get, put)
import           Data.Singletons         (SingI)
import qualified Data.Text               as T
import qualified Data.Vector             as V
import           Scientific.Workflow     hiding (_data)

import           Taiji.Core.Config       ()
import           Taiji.Types

type HiCWithSomeFile = HiC N [Either SomeFile (SomeFile, SomeFile)]

getHiCLoops :: [HiCWithSomeFile] -> [HiC S (File '[ChromosomeLoop] 'Bed)]
getHiCLoops inputs = concatMap split $ concatMap split $
    inputs & mapped.replicates.mapped.files %~ f
  where
    f fls = map fromSomeFile $
        filter (\x -> ChromosomeLoop `elem` getFileTags x) $ lefts fls

findActivePromoters :: SingI tags
                    => ATACSeq S (File tags 'NarrowPeak)
                    -> WorkflowConfig TaijiConfig (ATACSeq S (File tags 'Bed))
findActivePromoters input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks _taiji_output_dir >>= getPath
    let fun output fl = liftIO $ do
            peaks <- readBed' $ fl^.location :: IO [BED3]
            pro <- readPromoters anno
            writeBed' output $ findActivePromoters_ peaks pro
            return $ location .~ output $ emptyFile
    mapFileWithDefName (dir ++ "/") "_active_promoters.bed" fun input

findActivePromoters_ :: BEDLike b
                     => [b]        -- ^ Promoter activity indicator
                     -> [Promoter]
                     -> [Promoter]
findActivePromoters_ bed pro = runIdentity $ runConduit $ yieldMany pro .|
    intersectBed bed .| sinkList
{-# INLINE findActivePromoters_ #-}

-- | Get a list of potential TSS from GTF file
readPromoters :: FilePath -> IO [Promoter]
readPromoters = (fmap . concatMap) fn . readGenes
  where
    fn Gene{..} = map g $ nubSort tss
      where
        g x | geneStrand = BEDExt (asBed geneChrom (max 0 $ x - 5000) (x + 1000)) geneName
            | otherwise = BEDExt (asBed geneChrom (max 0 $ x - 1000) (x + 5000)) geneName
        tss | geneStrand = geneLeft : map fst geneTranscripts
            | otherwise = geneRight : map snd geneTranscripts
{-# INLINE readPromoters #-}

createLinkage :: ( ATACSeq S (File '[] 'Other)  -- ^ assignments
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ promoter activity
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ enhancer activity
                 )
              -> WorkflowConfig TaijiConfig (ATACSeq S (File '[] 'Other))
createLinkage (atac, pro, enh) = do
    dir <- asks ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir) >>= getPath
    let out = dir ++ "/linkages.bin"
    liftIO $ do
        activityPro <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_pro
        activityEnh <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_enh
        runResourceT $ runConduit $ readLinks fl_region activityPro activityEnh .|
            conduitPut put .| sinkFile out
    return $ atac & replicates.mapped.files .~ (emptyFile & location .~ out)
  where
    readLinks fl actP actE = sourceFileBS fl .| conduitGet2 get .| mapC f
      where
        f (nm, (ps, es)) =
            ( nm
            , M.toList $ fmap (mergeTFBS combineFn) $ M.fromListWith (++) $
                mapMaybe (g actP) ps
            , M.toList $ fmap (mergeTFBS combineFn) $ M.fromListWith (++) $
                mapMaybe (g actE) es) :: Linkage
        g act x = case getActivity x act of
            Nothing -> Nothing
            Just v -> Just ( getTFName x,
                [x & score.mapped %~ (transform_peak_height v *) . transform_site_pvalue] )
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
    fl_pro = either (^.location) (^.location) pro
    fl_enh = either (^.location) (^.location) enh
    fl_region = runIdentity (atac^.replicates) ^. files.location
    grp = atac^.groupName._Just

combineFn :: [Double] -> Double
combineFn = lp 2
{-# INLINE combineFn #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}

mergeTFBS :: ([Double] -> Double)  -- ^ Combining functin
          -> [BED] -> Double
mergeTFBS fn xs = fn $ runIdentity $ runConduit $
    mergeBedWith (maximum . map (fromJust . (^.score))) xs .| sinkList
{-# INLINE mergeTFBS #-}

-- | Assign activity to TFBS
getActivity :: BED -> BEDTree Double -> Maybe Double
getActivity tfbs beds = case IM.elems (intersecting beds tfbs) of
    [] -> Nothing
    xs -> Just $ maximum xs
{-# INLINE getActivity #-}

findTargets :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                         , File tag2 'Bed )        -- ^ TFBS
            -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
            -> WorkflowConfig TaijiConfig (ATACSeq S (File '[] 'Other))
findTargets atac hic = do
    dir <- asks ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir) >>= getPath
    let out = dir ++ "/TFBS_assignment.bin"
    liftIO $ do
        tfbs <- bedToTree (++) . map (\x -> (x, [x])) <$> readBed' (fl_tfbs^.location)
        promoters <- readBed' $ fl_pro^.location
        let enhancers2D = getRegulatoryDomains 1000000 promoters
        enhancers3D <- case hic of
            Nothing -> return []
            Just fl -> runResourceT $ runConduit $ read3DContact (fl^.location) .|
                getContactingRegions promoters .| sinkList
        let regions = groupRegions $
                zip (repeat Promoter) (mergeDomains promoters) ++
                zip (repeat Enhancer) (mergeDomains $ enhancers2D ++ enhancers3D)
        runResourceT $ runConduit $
            yieldMany regions .| findTargets_ tfbs .| conduitPut put .|
            sinkFile out
    return $ atac & replicates.mapped.files .~ (emptyFile & location .~ out)
  where
    [(fl_pro, fl_tfbs)] = atac^..replicates.folded.files
    grp = atac^.groupName._Just
    groupRegions = groupBy ((==) `on` ((^._data) . snd)) .
        sortBy (comparing ((^._data) . snd))

findTargets_ :: Monad m
             => BEDTree [BED]
             -> ConduitT [(DomainType, RegDomain)] (GeneName, ([BED], [BED])) m ()
findTargets_ tfbs = mapC ( \xs ->
    (snd (head xs) ^. _data, divide $ concat $ mapMaybe f xs) )
  where
    divide xs = let (a, b) = partition ((==Promoter) . fst) xs
                in (snd $ unzip a, snd $ unzip b)
    f (ty, region) = case IM.elems (intersecting tfbs region) of
        [] -> Nothing
        xs -> Just $ zip (repeat ty) $ concat xs
{-# INLINE findTargets_ #-}

-- | Merge 2D and 3D domains
mergeDomains :: [RegDomain] -> [RegDomain]
mergeDomains regions = runIdentity $ runConduit $ mergeBedWith f regions .|
    concatC .| sinkList
  where
    f xs = concatMap mergeFn $ groupBy ((==) `on` (^._data)) $
        sortBy (comparing (^._data)) xs
    mergeFn xs =
        let nm = head xs ^._data
            beds = runIdentity $ runConduit $ mergeBed xs .| sinkList
        in map (\x -> x & _data .~ nm) beds
{-# INLINE mergeDomains #-}

-- | Given a gene list , compute the rulatory domain for each gene
getRegulatoryDomains :: Int             -- ^ Extension length. A good default is 1M.
                     -> [Promoter] -- ^ A list of promoters
                     -> [RegDomain] -- ^ Regulatory domains
getRegulatoryDomains ext genes
    | null genes = error "No gene available for domain assignment!"
    | otherwise = concat $ loop $ [Nothing] ++ map Just basal ++ [Nothing]
  where
    loop (a:b:c:rest) = fn a b c : loop (b:c:rest)
    loop _            = []
    fn left (Just bed) right =
        [ bed & chromStart .~ leftPos & chromEnd .~ s
        , bed & chromStart .~ e & chromEnd .~ rightPos ]
      where
        chr = bed^.chrom
        s = bed^.chromStart
        e = bed^.chromEnd
        leftPos
            | isNothing left || chr /= fromJust left ^. chrom = max (s - ext) 0
            | otherwise = min s $ max (s - ext) $ fromJust left ^. chromEnd
        rightPos
            | isNothing right || chr /= fromJust right ^. chrom = e + ext   -- TODO: bound check
            | otherwise = max e $ min (e + ext) $ fromJust right ^. chromStart
    fn _ _ _ = undefined
    basal = V.toList $ fromSorted $ sortBed genes
{-# INLINE getRegulatoryDomains #-}

-- | Read 3D contacts from a file, where each line contains 6 fields separated
-- by Tabs corresponding to 2 interacting loci. Example:
-- chr1 [TAB] 11 [TAB] 100 [TAB] chr2 [TAB] 23 [TAB] 200
read3DContact :: FilePath -> ConduitT i (BED3, BED3) (ResourceT IO) ()
read3DContact input = sourceFileBS input .| linesUnboundedAsciiC .| mapC f .|
    filterC g
  where
    f x = let (chr1:s1:e1:chr2:s2:e2:_) = B.split '\t' x
          in (asBed chr1 (readInt s1) (readInt e1), asBed chr2 (readInt s2) (readInt e2))
    g (x1,x2) | x1^.chrom /= x2^.chrom = True
              | x1^.chromEnd < x2^.chromStart = x2^.chromStart - x1^.chromEnd > 5000
              | x2^.chromEnd < x1^.chromStart = x1^.chromStart - x2^.chromEnd > 5000
              | otherwise = False
{-# INLINE read3DContact #-}

-- | Retrieve regions that contact with the gene's promoter.
getContactingRegions :: Monad m
                     => [Promoter]
                     -> ConduitT (BED3, BED3) RegDomain m ()
getContactingRegions promoters =
    let basal = bedToTree (++) $ flip map promoters $
            \pro -> (pro^._bed, [pro^._data])
     in concatMapC $ \(locA, locB) -> map (BEDExt locB) (intersect basal locA) ++
            map (BEDExt locA) (intersect basal locB)
  where
    intersect t x = nubSort $ concat $ IM.elems $ intersecting t x
{-# INLINE getContactingRegions #-}
