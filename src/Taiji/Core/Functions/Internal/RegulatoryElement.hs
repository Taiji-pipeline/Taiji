{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
module Taiji.Core.Functions.Internal.RegulatoryElement where

import           Bio.Data.Bed
import           Bio.RealWorld.GENCODE
import           Bio.Utils.Misc          (readInt)
import           Conduit
import           Control.Lens            ((^.), (.~), (&))
import qualified Data.ByteString.Char8   as B
import           Data.CaseInsensitive    (CI)
import qualified Data.HashMap.Strict     as M
import qualified Data.IntervalMap.Strict as IM
import           Data.List.Ordered       (nubSort)
import           Data.Maybe              (fromJust, isNothing, mapMaybe)
import qualified Data.Vector             as V

import           Taiji.Types

createLinkage :: [TFBS] -> [(CI B.ByteString, [(CI B.ByteString, [TFBS])])]
createLinkage tfbs = M.toList $ fmap
    (M.toList . M.fromListWith (++) . map (\x -> (_tf_name $ _ext_data x, [x]))) $
    M.fromListWith (++) $ mapMaybe f tfbs
  where
    f x = let gene = _target_gene $ _ext_data x
          in case gene of
              Nothing -> Nothing
              Just g  -> Just (g, [x])

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

findActivePromoters :: BEDLike b
                    => [b]        -- ^ Promoter activity indicator
                    -> [Promoter]
                    -> [Promoter]
findActivePromoters bed pro = runIdentity $ runConduit $ yieldMany pro .|
    intersectBed bed .| sinkList

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

-- | Given a gene list and the rule, compute the rulatory domain for each gene
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
