{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser   (readHiC, readHiCTSV)
import           Control.Lens
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Reader         (asks)
import qualified Data.Map.Strict              as M
import           Data.Maybe                   (fromJust)
import           Data.Monoid                  ((<>))
import           Scientific.Workflow

import           Taiji.Core.Network
import           Taiji.Core.Ranking
import           Taiji.Core.RegulatoryElement
import           Taiji.Types                  (_taiji_input)

builder :: Builder ()
builder = do
    nodePS 1 "Find_Active_Promoter" 'findActivePromoters $ do
        note .= "Identify active promoters. Promoters are defined by " <>
            "-5000 ~ +1000 regions around annotated transcription start sites. " <>
            "If a promoter is overlapped with ATAC-seq peaks, we assume it is active."

    nodeS "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- if ".tsv" == reverse (take 4 $ reverse input)
                then readHiCTSV input "HiC"
                else readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ do
            submitToRemote .= Just False
            note .= "Read HiC loops from input file."

    node' "Create_Linkage_Prep" [| \(activePro, tfbs, atac_peaks, chip_peaks, hic, expr) ->
        let getFile x = (x^.groupName._Just, x^.replicates._2.files)
            tfbsMap = M.fromList $ map getFile tfbs
            atacFileMap = fmap Left $ M.fromList $ map getFile atac_peaks
            chipFileMap = M.fromList $ map getFile chip_peaks
            hicFileMap = M.fromList $ map getFile hic
        in flip map activePro $ \e ->
            let grp = e^.groupName._Just
                pro = M.findWithDefault undefined grp atacFileMap
                enh = M.findWithDefault pro grp chipFileMap
                hic' = M.lookup grp hicFileMap
                e' = e & replicates.mapped.files %~ (\f -> (f, fromJust $ M.lookup grp tfbsMap))
            in (e', pro, enh, hic', expr)
        |] $ do
            note .= "Prepare for parallel execution."
            submitToRemote .= Just False
    nodePS 1 "Create_Linkage" 'createLinkage $ do
        remoteParam .= "--mem=20000 -p gpu"
    path ["Create_Linkage_Prep", "Create_Linkage"]

    node' "Compute_Ranks_Prep" [| \(x, expr) -> zip x $ repeat expr |] $ do
        submitToRemote .= Just False
    nodePS 1 "Compute_Ranks" 'computeRanks $ do
        note .= "Perform personalized Pagerank."
        remoteParam .= "--mem=20000 -p gpu"
    nodeS "Output_Ranks" 'outputRanks $ return ()

    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Ranks"]

