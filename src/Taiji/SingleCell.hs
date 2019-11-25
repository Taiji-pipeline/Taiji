{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.SingleCell ( builder ) where

import Control.Workflow

import Taiji.SingleCell.Cluster

builder :: Builder ()
builder = do
    --node "Compute_Ranks_SC_Prep" 'prepDataSet $ return ()
    {-
    nodePS 1 "Compute_Ranks_SC" 'computeRanksSC $ return ()
    path ["Compute_Ranks_SC_Prep", "Compute_Ranks_SC"]
    -}

    node "Compute_Ranks_SC_Prep" [| return . prepareData |] $ return ()
    nodePar "Compute_Ranks_SC" 'computeRanksCluster $ return ()
    node "Output_Ranks_SC" [| outputRanksCluster "/Rank/Cluster/" |] $ return ()
    path ["Compute_Ranks_SC_Prep", "Compute_Ranks_SC", "Output_Ranks_SC"]