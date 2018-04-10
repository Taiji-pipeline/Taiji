{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser   (readHiC)
import           Control.Lens
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Reader         (asks)
import qualified Data.Map.Strict              as M
import           Data.Maybe                   (fromJust)
import           Data.Monoid                  ((<>))
import           Scientific.Workflow

import           Taiji.Core.Network
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
            hic <- readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ do
            submitToRemote .= Just False
            note .= "Read HiC loops from input file."

    node' "Find_TF_Target_Prep" [| \(activePro, tfbs, peaks, hic) ->
        let mkDict xs = M.fromList $ map (\x ->
                (x^.groupName, head $ x^..replicates.folded.files)) xs
            lookup' x xs = M.lookup (x^.groupName) xs
            tfbs' = mkDict tfbs
            hic' = mkDict hic
            peaks' = mkDict peaks
        in flip map activePro $ \e ->
            ( e & replicates.mapped.files %~ (\f -> (f, fromJust $ lookup' e tfbs'))
            , fromJust $ lookup' e peaks'
            , fromJust $ lookup' e peaks'
            , M.lookup (e^.groupName) hic' )
        |] $ note .= "Prepare for parallel execution."

    nodePS 1 "Find_TF_Target" [| \(x1,x2,x3,x4) -> findTargets x1 x2 x3 x4 |] $ do
        note .= "Assign TFs to their target genes. We borrow the concept of " <>
                "gene regulatory domain from GREAT. Gene regulatory domain " <>
                "definition: Active promoters are used to define the basal " <>
                "regulatory domains of genes. The gene regulatory domain is " <>
                "extended in both directions to the nearest gene's basal domain " <>
                "but no more than the maximum extension in one direction." <>
                "TF binding sites located in gene regulatory domains are then " <>
                "assigned to corresponding genes."
    path ["Find_TF_Target_Prep", "Find_TF_Target"]

    node' "Compute_Ranks_Prep" [| \(links, expr) -> zip links $ repeat expr
        |] $ note .= "Prepare for parallel execution."
    nodePS 1 "Compute_Ranks" 'computeRanks $ do
        note .= "Perform personalized Pagerank."
        remoteParam .= "--mem=20000 -p gpu"
    nodeS "Output_Rank" 'outputRank $ return ()
    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Rank"]
