{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser (readHiC)
import           Control.Lens
import qualified Data.Map.Strict      as M
import           Data.Monoid          ((<>))
import           Scientific.Workflow
import           Control.Monad.IO.Class                (liftIO)
import           Control.Monad.Reader                  (asks)

import           Taiji.Core.Exporter  (exportResults)
import           Taiji.Core.Functions
import           Taiji.Types (_taiji_input)

builder :: Builder ()
builder = do
    nodePS 1 "Find_Active_Promoter" 'getActivePromoter $ do
        note .= "Identify active promoters. Promoters are defined by " <>
            "-5000 ~ +1000 regions around annotated transcription start sites. " <>
            "If a promoter is overlapped with ATAC-seq peaks, we assume it is active."

    nodeS "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ do
            submitToRemote .= Just False
            note .= "Read HiC loops from input file."

    node' "Link_Gene_TF_Prep" [| \(activePro, tfbs, hic) ->
        let getFile e = head $ e^..replicates.folded.files
            tfbs' = M.fromList $ map (\x -> (x^.groupName, getFile x)) tfbs
            hic' = M.fromList $ map (\x -> (x^.groupName, x)) hic
        in flip map activePro $ \e ->
             let a = M.findWithDefault undefined (e^.groupName) tfbs'
             in ( e & replicates.mapped.files %~ (\f -> (f,a))
                , M.lookup (e^.groupName) hic' )
        |] $ note .= "Prepare for parallel execution."
    nodePS 1 "Link_Gene_TF" [| \(x, y) -> linkGeneToTFs x y |] $ do
        note .= "Assign TFs to their target genes. We borrow the concept of " <>
            "gene regulatory domain from GREAT. Gene regulatory domain " <>
            "definition: Active promoters are used to define the basal " <>
            "regulatory domains of genes. The gene regulatory domain is " <>
            "extended in both directions to the nearest gene's basal domain " <>
            "but no more than the maximum extension in one direction." <>
            "TF binding sites located in gene regulatory domains are then " <>
            "assigned to corresponding genes."
    path ["Link_Gene_TF_Prep", "Link_Gene_TF"]

    node' "TFRank_Prep" [| \(a,b) -> zip (repeat a) b |] $ note .= "Prepare for parallel execution."
    nodePS 1 "TFRank" 'getTFRanks $ do
        note .= ""
    nodeS "Output_Rank" 'outputRank $ return ()
    path ["TFRank_Prep", "TFRank", "Output_Rank"]

    nodeS "Export_Results" 'exportResults $ return ()
