{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.RegulatoryElement
    ( getHiCLoops
    , findActivePromoters
    , findTargets
    ) where

import           Bio.Data.Bed
import           Bio.Data.Bed.Types (BED3(..), fromSorted)
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils      (getPath)
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc          (readInt)
import           Conduit
import Control.Arrow (second)
import           Control.Lens
import           Control.Monad.Reader    (asks)
import qualified Data.ByteString.Char8   as B
import           Data.Either             (lefts)
import           Data.Function           (on)
import qualified Data.IntervalMap.Generic.Strict as IM
import qualified Data.HashMap.Strict as M
import           Data.List               
import           Data.List.Ordered       (nubSort)
import           Data.Maybe              (fromJust, isNothing)
import           Data.Monoid             ((<>))
import           Data.Ord
import           Data.Singletons         (SingI)
import qualified Data.Text               as T
import qualified Data.Vector             as V
import           Scientific.Workflow     hiding (_data)
import           Text.Printf             (printf)

import           Taiji.Types
import           Taiji.Core.Utils

getHiCLoops :: [HiC N [Either SomeFile (SomeFile, SomeFile)]]
            -> [HiC S (File '[ChromosomeLoop] 'Bed)]
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
    dir <- asks ((<> "/Promoters") . _taiji_output_dir) >>= getPath
    let output = printf "%s/%s_rep%d_active_promoters.bed" dir
            (T.unpack $ input^.eid) (input^.replicates._1)
    input & replicates.traverse.files %%~ liftIO . ( \fl -> do
        peaks <- readBed' $ fl^.location :: IO [BED3]
        pro <- readPromoters anno
        writeBed' output $ findActivePromoters_ peaks pro
        return $ location .~ output $ emptyFile
        )

findActivePromoters_ :: BEDLike b
                     => [b]        -- ^ Promoter activity indicator
                     -> [Promoter]
                     -> [Promoter]
findActivePromoters_ bed pro = runIdentity $ runConduit $ yieldMany pro .|
    intersectBed bed .| sinkList
{-# INLINE findActivePromoters_ #-}


findTargets :: MonadResource m
            => BEDTree [SiteInfo]   -- ^ TF binding sites
            -> File tag1 'Bed         -- ^ Active promoters
            -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
            -> ConduitT () (GeneName, ([TFBS], [TFBS])) m ()
findTargets tfbs fl_pro hic = do
    promoters <- liftIO $ readBed' $ fl_pro^.location
    let chromLoops = case hic of
            Nothing -> return ()
            Just fl -> read3DContact (fl^.location)
    chromLoops .| findTargets_ tfbs promoters


findTargets_ :: Monad m
             => BEDTree [SiteInfo]   -- ^ TF binding sites
             -> [Promoter] -- ^ A list of promoters
             -> ConduitT (BED3, BED3)  -- Chromatin loops
                    (GeneName, ([TFBS], [TFBS])) m ()
findTargets_ tfbs promoters = getRegulaDomain promoters .|
    mapC (second (divide . concatMap f))
  where
    divide xs = let (a, b) = partition ((==Promoter) . fst) xs
                in (snd $ unzip a, snd $ unzip b)
    f region = zip (repeat $ region^._data) $ concatMap mkTFBS $ IM.toList $
        intersecting tfbs $ region^._bed
      where
        mkTFBS (x, s) = let bed = BED3 (region^.chrom) (IM.lowerBound x)
                                (IM.upperBound x)
                        in zipWith BEDExt (repeat bed) s
{-# INLINE findTargets_ #-}

-- | The regulatory domain of a gene is the genomic regions possessing
-- regulatory effects on the gene, such as promoters and enhancers.
-- To find the regulatory domain, distal enhancers are first identified
-- according to chromatin interactions. Then the promoter 
getRegulaDomain :: Monad m
                => [Promoter] -- ^ A list of promoters
                -> ConduitT
                       (BED3, BED3)  -- Chromatin loops
                       (GeneName, [BEDExt BED3 DomainType])   -- Regulatory domains for each gene
                       m ()
getRegulaDomain inputRegions = do
    enhancers3D <- getDistalDomain3D inputRegions .| sinkList
    let enhancers = map (\x -> (x^._data, [BEDExt (x^._bed) Enhancer])) $
            mergeDomains $ enhancers2D ++ enhancers3D
        promoters = map (\x -> (x^._data, [BEDExt (x^._bed) Promoter])) $
            mergeDomains inputRegions
    yieldMany $ M.toList $ M.fromListWith (++) $ enhancers ++ promoters
  where
    enhancers2D = getDistalDomain2D 1000000 inputRegions
{-# INLINE getRegulaDomain #-}


-- | Merge 2D and 3D domains
mergeDomains :: [BEDExt BED3 GeneName] -> [BEDExt BED3 GeneName]
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

-- | Compute the 2D distal domains (enhancer regions) for a list of genes.
-- The regulatory domain doesn't cover genes' promoters.
getDistalDomain2D :: Int             -- ^ Extension length. A good default is 1M.
                  -> [Promoter] -- ^ A list of promoters
                  -> [RegDomain] -- ^ Regulatory domains
getDistalDomain2D ext genes
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
{-# INLINE getDistalDomain2D #-}

-- | Compute the 3D distal domains (enhancer regions) for a list of genes.
-- The regulatory domain doesn't cover genes' promoters.
-- Input: Chromatin loops represented by pairs of interactng loci.
getDistalDomain3D :: Monad m
                  => [Promoter]   -- ^ Gene promoters
                  -> ConduitT (BED3, BED3) RegDomain m ()
getDistalDomain3D promoters =
    let basal = bedToTree (++) $ flip map promoters $
            \pro -> (pro^._bed, [pro^._data])
     in concatMapC $ \(locA, locB) -> map (BEDExt locB) (intersect basal locA) ++
            map (BEDExt locA) (intersect basal locB)
  where
    intersect t x = nubSort $ concat $ IM.elems $ intersecting t x
{-# INLINE getDistalDomain3D #-}


-------------------------------------------------------------------------------
-- IO
-------------------------------------------------------------------------------

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

-- | Read 3D contacts from a file, where each line contains 6 fields separated
-- by Tabs corresponding to 2 interacting loci. Example:
-- chr1 [TAB] 11 [TAB] 100 [TAB] chr2 [TAB] 23 [TAB] 200
read3DContact :: MonadResource m
              => FilePath -> ConduitT i (BED3, BED3) m ()
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

