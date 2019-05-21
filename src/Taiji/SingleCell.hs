{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.SingleCell ( builder ) where

import Taiji.Prelude
import Taiji.SingleCell.Cell
import Taiji.SingleCell.Cluster

builder :: Builder ()
builder = do
    node "Compute_Ranks_SC_Prep" 'prepDataSet $ submitToRemote .= Just False
    nodePS 1 "Compute_Ranks_SC" 'computeRanksSC $ return ()
    path ["Compute_Ranks_SC_Prep", "Compute_Ranks_SC"]

    node' "Compute_Ranks_SC_Cluster_Prep" [| \(x, y) -> zip (repeat x) y |] $ return ()
    nodePS 1 "Compute_Ranks_SC_Cluster" 'computeRanksCluster $ return ()
    path ["Compute_Ranks_SC_Cluster_Prep", "Compute_Ranks_SC_Cluster"]
