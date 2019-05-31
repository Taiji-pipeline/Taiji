{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.SingleCell ( builder ) where

import Control.Workflow

import Taiji.SingleCell.Cell
import Taiji.SingleCell.Cluster

builder :: Builder ()
builder = do
    node "Compute_Ranks_SC_Prep" 'prepDataSet $ return ()
    {-
    nodePS 1 "Compute_Ranks_SC" 'computeRanksSC $ return ()
    path ["Compute_Ranks_SC_Prep", "Compute_Ranks_SC"]
    -}

    node "Compute_Ranks_SC_Cluster_Prep" [| return . prepareData |] $ return ()
    nodePar "Compute_Ranks_SC_Cluster" 'computeRanksCluster $ return ()
    node "Output_Ranks_SC_Cluster" 'outputRanksCluster $ return ()
    path [ "Compute_Ranks_SC_Cluster_Prep", "Compute_Ranks_SC_Cluster",
        "Output_Ranks_SC_Cluster" ]
