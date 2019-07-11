{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Core (builder) where

import           Bio.Pipeline (getPath)
import           Bio.Data.Experiment.Parser   (readHiC, readHiCTSV)
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Map.Strict              as M
import Control.Workflow

import           Taiji.Core.Network
import           Taiji.Core.Ranking
import           Taiji.Core.RegulatoryElement
import Taiji.Prelude

aggregate :: MonadIO m
          => ( [ATACSeq S f1]   -- ^ TFBS
             , [ATACSeq S f2]   -- ^ peaks for promoter
             , [HiC S f3]  -- ^ HiC loops
             , Maybe (File '[] 'Tsv)       -- ^ Expression
             , Maybe (File '[] 'Tsv) )         -- ^ estimated Expression from ATAC-seq
          -> m [ (ATACSeq S f1, f2, Maybe f3, Maybe (File '[] 'Tsv)) ]
aggregate (tfbs, atac_peaks, hic, rnaE, atacE) = liftIO $ do
    grps <- case expr of
        Nothing -> return Nothing
        Just fl -> Just . map (T.pack . B.unpack) . tail . B.split '\t' .
            head . B.lines <$> B.readFile (fl^.location)
    return $ flip mapMaybe tfbs $ \e ->
        let grp = e^.groupName._Just
            pro = M.findWithDefault undefined grp atacFileMap
            hic' = M.lookup grp hicFileMap
        in case grps of
            Nothing -> Just (e, pro, hic', expr)
            Just grps' -> if grp `elem` grps'
                then Just (e, pro, hic', expr)
                else Nothing
  where
    expr = maybe atacE Just rnaE
    atacFileMap = M.fromList $ map getFile atac_peaks
    hicFileMap = M.fromList $ map getFile hic
    getFile x = (x^.groupName._Just, x^.replicates._2.files)

builder :: Builder ()
builder = do
    node "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- if ".tsv" == reverse (take 4 $ reverse input)
                then readHiCTSV input "HiC"
                else readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ doc .= "Read HiC loops from input file."

    node "Create_Linkage_Prep" 'aggregate $
        doc .= "Prepare for parallel execution."
    nodePar "Create_Linkage" 'saveAssociations $ memory .= 20
    path ["Create_Linkage_Prep", "Create_Linkage"]

    node "Compute_Ranks_Prep" [| \(x, expr) -> return $ zip x $ repeat expr |] $
        return ()
    nodePar "Compute_Ranks" 'computeRanks $ do
        doc .= "Perform personalized Pagerank."
        memory .= 20
    node "Output_Ranks" [| \input -> do
        dir <- asks _taiji_output_dir >>= getPath
        let output1 = dir ++ "/GeneRanks.tsv"
            output2 = dir ++ "/GeneRanks_PValues.tsv"
            output3 = dir ++ "/GeneRanks.html"
        liftIO $ outputRanks output1 output2 output3 input
        |] $ return ()
    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Ranks"]