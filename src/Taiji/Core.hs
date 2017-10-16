{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Control.Lens
import qualified Data.Map.Strict      as M
import           Scientific.Workflow

import           Taiji.Core.Functions
import           Taiji.Core.Exporter (exportResults)

builder :: Builder ()
builder = do
    nodePS 1 "Find_Active_Promoter" 'getActivePromoter $ return ()

    node' "Link_Gene_TF_Prep" [| \(activePro, tfbs) ->
        let getFile e = head $ e^..replicates.folded.files
            tfbs' = M.fromList $ map (\x -> (x^.groupName, getFile x)) tfbs
        in flip map activePro $ \e ->
             let a = M.findWithDefault undefined (e^.groupName) tfbs'
             in e & replicates.mapped.files %~ (\f -> (f,a))
        |] $ return ()
    nodePS 1 "Link_Gene_TF" 'linkGeneToTFs $ return ()
    path ["Link_Gene_TF_Prep", "Link_Gene_TF"]

    node' "TFRank_Prep" [| \(a,b) -> zip (repeat a) b |] $ return ()
    nodeP 1 "TFRank" 'getTFRanks $ return ()
    nodeS "Output_Rank" 'outputRank $ return ()
    path ["TFRank_Prep", "TFRank", "Output_Rank"]

    nodeS "Export_Results" 'exportResults $ return ()
