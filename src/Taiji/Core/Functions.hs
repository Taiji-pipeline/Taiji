{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Core.Functions
    ( getActivePromoter
    ) where

import           Bio.Data.Bed           (BED (..), BED3, BEDLike (..),
                                         intersectBed, readBed', writeBed')
import           Bio.Data.Experiment
import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.NGS
import           Bio.Pipeline.Utils     (getPath)
import           Bio.Utils.Misc         (readInt)
import           Conduit
import           Control.Lens
import           Control.Monad.Reader   (asks)
import qualified Data.ByteString.Char8  as B
import           Data.List.Ordered      (nubSort)
import           Data.Maybe             (fromJust)
import           Data.Singletons        (SingI)
import           Scientific.Workflow
import           System.IO.Temp         (withTempFile)

import           Taiji.Core.Config

getActivePromoter :: SingI tags
                  => ATACSeq (File tags 'Bed)
                  -> WorkflowConfig TaijiConfig (ATACSeq (File tags 'Bed))
getActivePromoter input = do
    anno <- fromJust <$> asks _taiji_annotation
    dir <- asks _taiji_output_dir >>= getPath
    let fun output fl = liftIO $ withTempFile "./" "tmp_macs2_file." $ \tmp _ -> do
            _ <- callPeaks tmp fl Nothing $ cutoff .= QValue 0.1
            peaks <- readBed' tmp :: IO [BED3]
            tss <- getActiveTSS anno peaks
            writeBed' output tss
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
        runIdentity $ yieldMany promoters =$= intersectBed peaks $$ sinkList
  where
    f xs = let name = B.filter (/='\"') $ last $ B.split ' ' $ head $
                      filter (B.isInfixOf "gene_name") $ B.split ';' $ last xs
               start = readInt $ xs !! 3
               end = readInt $ xs !! 4
               strand = if xs !! 6 == "+" then True else False
          in ((xs !! 0, if strand then start else end, strand), name)
    g xs = not $ B.isPrefixOf "#" (head xs) || (xs !! 2 /= "transcript")
{-# INLINE getActiveTSS #-}
