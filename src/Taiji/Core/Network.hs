{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network
    ( combineNetwork
    , mkNetwork
    , readNodesAndEdges
    , createLinkage
    , getExternalLinks
    ) where

import Data.List.Ordered (nubSortBy)
import Data.Ord (comparing)
import qualified Data.HashMap.Strict as M
import IGraph

import           Taiji.Core.Config                 ()
import           Taiji.Types
import Taiji.Core.Network.DeNovo
import Taiji.Core.Network.External

combineNetwork :: ([NetNode], [((GeneName, GeneName), Double)])
               -> ([NetNode], [((GeneName, GeneName), Double)])
               -> Graph 'D NetNode Double
combineNetwork (nds1, es1) (nds2, es2) = mkGraph nds es
  where
    nds = nubSortBy (comparing _node_name) $ nds1 ++ nds2
    nodeMap = M.fromList $ zip (map _node_name nds) [0..]
    lookupId x = M.lookupDefault undefined x nodeMap
    es = M.toList $ M.fromListWith max $
        map (\((a, b), x) -> ((lookupId b, lookupId a), x)) $ es1 ++ es2
