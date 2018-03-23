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
import           Bio.Data.Experiment
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils       (asDir, getPath)
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc           (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.Reader     (asks)
import           Control.Monad.State.Lazy (StateT, execStateT, modify)
import qualified Data.ByteString.Char8    as B
import           Data.CaseInsensitive     (mk)
import           Data.Either              (lefts)
import qualified Data.IntervalMap.Strict  as IM
import           Data.List.Ordered        (nubSort)
import           Data.Maybe               (fromJust, isNothing)
import           Data.Singletons          (SingI)
import qualified Data.Vector              as V
import           Scientific.Workflow      hiding (_data)

import           Taiji.Core.Config        ()
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
    dir <- asks (asDir . _taiji_output_dir) >>= getPath
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


findTargets :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                         , File tag2 'Bed )        -- ^ TFBS
            -> File tag3 'NarrowPeak
            -> File tag4 'NarrowPeak
            -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
            -> IO [TFBS]
findTargets atac activityPro activityEnh hic = do
    promoters <- readBed' $ fl_pro^.location
    let enhancers = getRegulatoryDomains 1000000 promoters
    distalEhancers <- case hic of
        Nothing -> return []
        Just fl -> runResourceT $ runConduit $ read3DContact (fl^.location) .|
            getContactingRegions promoters .| sinkList
    activityPro' <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
        readBed' $ activityPro^.location
    activityEnh' <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
        readBed' $ activityEnh^.location

    let action = runConduit $ readBed (fl_tfbs^.location) .|
            findTargets_ promoters enhancers distalEhancers activityPro' activityEnh' .|
            sinkNull
    execStateT action []
  where
    [(fl_pro, fl_tfbs)] = atac^..replicates.folded.files

findTargets_ :: Monad m
             => [Promoter]
             -> [RegDomain]
             -> [RegDomain]
             -> BEDTree Double   -- ^ promoter activity
             -> BEDTree Double   -- ^ enhancer activity
             -> ConduitT BED TFBS (StateT [TFBS] m) ()
findTargets_ promoters enhancers distalEhancers activityPro activityEnh= mapC toSite .|
    assignTFBS promoters .| concatMapMC (fn Pro activityPro) .|
    assignTFBS distalEhancers .| concatMapMC (fn NearbyEnh activityEnh) .|
    assignTFBS enhancers .| concatMapMC (fn DistalEnh activityEnh)
  where
    toSite x = BEDExt (x :: BED) $ SiteInfo (mk $ fromJust $ x^.name)
        Nothing Nothing Nothing
    fn _ _ (Left x) = return $ Just x
    fn t tree (Right x) = case getActivity x tree of
        Nothing -> return Nothing
        Just a -> do
            let x' = x & _data.target_through .~ Just t
                       & _data.peak_signal .~ Just a
            modify (x':)
            return Nothing
{-# INLINE findTargets_ #-}

assignTFBS :: Monad m
           => [RegDomain]   -- ^ Regulatory element
           -> ConduitT TFBS (Either TFBS TFBS) m ()
assignTFBS regions = intersectBedWith f regions .| concatC
  where
    f tfbs promoters
        | null promoters = [Left tfbs]
        | otherwise = zipWith (\a b -> Right $ a &
            _data.target_gene .~ Just (b^._data))
            (repeat tfbs) promoters
{-# INLINE assignTFBS #-}

-- | Assign activity to TFBS
getActivity :: TFBS -> BEDTree Double -> Maybe Double
getActivity tfbs beds = case IM.elems (intersecting beds tfbs) of
    [] -> Nothing
    xs -> Just $ maximum xs
{-# INLINE getActivity #-}

-- | Given a gene list , compute the rulatory domain for each gene
getRegulatoryDomains :: Int             -- ^ Extension length. A good default is 1M.
                     -> [Promoter] -- ^ A list of promoters
                     -> [RegDomain] -- ^ Regulatory domains
getRegulatoryDomains ext genes
    | null genes = error "No gene available for domain assignment!"
    | otherwise = loop $ [Nothing] ++ map Just basal ++ [Nothing]
  where
    loop (a:b:c:rest) = fn a b c : loop (b:c:rest)
    loop _            = []
    fn left (Just bed) right = bed & chromStart .~ leftPos & chromEnd .~ rightPos
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
read3DContact input = sourceFileBS input .| linesUnboundedAsciiC .| mapC f
  where
    f x = let (chr1:s1:e1:chr2:s2:e2:_) = B.split '\t' x
          in (asBed chr1 (readInt s1) (readInt e1), asBed chr2 (readInt s2) (readInt e2))
{-# INLINE read3DContact #-}

-- | Retrieve regions that contact with the gene's promoter.
getContactingRegions :: Monad m
                     => [Promoter]
                     -> ConduitT (BED3, BED3) RegDomain m ()
getContactingRegions promoters =
    let basal = bedToTree (++) $ flip map promoters $
            \pro -> (_ext_bed pro, [_ext_data pro])
     in concatMapC $ \(locA, locB) -> map (BEDExt locB) (intersect basal locA) ++
            map (BEDExt locA) (intersect basal locB)
  where
    intersect t x = nubSort $ concat $ IM.elems $ intersecting t x
{-# INLINE getContactingRegions #-}
