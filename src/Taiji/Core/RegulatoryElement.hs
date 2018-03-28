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
import           Bio.Pipeline.Utils       (asDir, getPath)
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc           (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.Reader     (asks)
import           Control.Monad.State.Lazy (StateT, execStateT, get, put)
import qualified Data.ByteString.Char8    as B
import           Data.CaseInsensitive     (CI, mk)
import           Data.Either              (lefts)
import qualified Data.HashMap.Strict      as M
import qualified Data.IntervalMap.Strict  as IM
import           Data.List.Ordered        (nubSort)
import           Data.Maybe               (fromJust, fromMaybe, isNothing)
import           Data.Singletons          (SingI)
import qualified Data.Vector              as V
import           Scientific.Workflow      hiding (_data)

import           Taiji.Core.Config        ()
import           Taiji.Types

-- | Gene and its regulators
type GeneName = CI B.ByteString
type Linkage = (GeneName, [(GeneName, (Maybe (Double, Int), Maybe (Double, Int), Maybe (Double, Int)))])

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

createLinkage :: St -> [Linkage]
createLinkage (pro, enh, distal) = M.toList $ combineHashMaps fn pro enh distal
  where
    fn a b c = M.toList $ combineHashMaps g (fromMaybe M.empty a) (fromMaybe M.empty b) (fromMaybe M.empty c)
    g a b c = (a, b, c)

    f x = let gene = _target_gene $ _ext_data x
          in case gene of
              Nothing -> Nothing
              Just g  -> Just (g, [x])
{-# INLINE createLinkage #-}

findTargets :: ATACSeq S ( File tag1 'Bed          -- ^ Active promoters
                         , File tag2 'Bed )        -- ^ TFBS
            -> File tag3 'NarrowPeak
            -> File tag4 'NarrowPeak
            -> Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
            -> IO St
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
    execStateT action (M.empty, M.empty, M.empty)
  where
    [(fl_pro, fl_tfbs)] = atac^..replicates.folded.files

type St = ( M.HashMap GeneName (M.HashMap GeneName (Double, Int))
          , M.HashMap GeneName (M.HashMap GeneName (Double, Int))
          , M.HashMap GeneName (M.HashMap GeneName (Double, Int)) )

combineHashMaps :: (Maybe v -> Maybe v -> Maybe v -> v')
                -> M.HashMap GeneName v
                -> M.HashMap GeneName v
                -> M.HashMap GeneName v
                -> M.HashMap GeneName v'
combineHashMaps fn m1 m2 m3 = M.fromList $ flip map keys $
    \k -> (k, fn (M.lookup k m1) (M.lookup k m2) (M.lookup k m3))
  where
    keys = nubSort $ M.keys m1 ++ M.keys m2 ++ M.keys m3

combineFunc :: Double -> Double -> Double
combineFunc = max

findTargets_ :: Monad m
             => [Promoter]
             -> [RegDomain]
             -> [RegDomain]
             -> BEDTree Double   -- ^ promoter activity
             -> BEDTree Double   -- ^ enhancer activity
             -> ConduitT BED BED (StateT St m) ()
findTargets_ promoters enhancers distalEhancers activityPro activityEnh=
    assignTFBS promoters .| concatMapMC (fn 1 activityPro) .|
    assignTFBS distalEhancers .| concatMapMC (fn 2 activityEnh) .|
    assignTFBS enhancers .| concatMapMC (fn (3::Int) activityEnh)
  where
    fn _ _ (Left x) = return $ Just x
    fn idx tree (Right (tf, gene)) = case getActivity tf tree of
        Nothing -> return Nothing
        Just v -> do
            (a,b,c) <- get
            let tfName = mk $ tf^.name._Just
                sc = transform_site_pvalue (fromJust $ tf^.score) *
                    transform_peak_height v
                table = case idx of
                    1 -> a
                    2 -> b
                    3 -> c
                    _ -> undefined
                table' =
                    let f (Just m) = Just $ M.alter g tfName m
                        f Nothing  = Just $ M.singleton tfName (sc, 1)
                        g Nothing         = Just (sc, 1)
                        g (Just (x1, x2)) = Just (combineFunc x1 sc, x2 + 1)
                    in M.alter f gene table
                newSt = case idx of
                    1 -> (table', b, c)
                    2 -> (a, table', c)
                    3 -> (a, b, table')
                    _ -> undefined
            put newSt
            return Nothing
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE findTargets_ #-}

assignTFBS :: Monad m
           => [RegDomain]   -- ^ Regulatory element
           -> ConduitT BED (Either BED (BED, CI B.ByteString)) m ()
assignTFBS regions = intersectBedWith f regions .| concatC
  where
    f tfbs promoters
        | null promoters = [Left tfbs]
        | otherwise = zipWith (\a b -> Right (a,  b^._data)) (repeat tfbs) promoters
{-# INLINE assignTFBS #-}

-- | Assign activity to TFBS
getActivity :: BED -> BEDTree Double -> Maybe Double
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
