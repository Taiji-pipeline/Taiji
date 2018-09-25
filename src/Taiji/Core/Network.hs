{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network
    ( combineNetwork
    , mkNetwork
    , readNetwork
    , createLinkage
    , getExternalLinks
    ) where

import Control.Arrow ((&&&))
import Data.List.Ordered (nubSortBy)
import Data.Ord (comparing)
import qualified Data.HashMap.Strict as M
import qualified Data.ByteString.Char8 as B
import IGraph
import           Bio.Utils.Misc                    (readDouble)
import Conduit
import           Data.CaseInsensitive              (mk)

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

-- | Build the network from files containing the information of nodes and edges.
mkNetwork :: FilePath  -- ^ nodes
          -> FilePath  -- ^ edges
          -> IO (Graph 'D NetNode Double)
mkNetwork nodeFl edgeFl = do
    nodeMap <- M.fromList . map ((_node_name &&& id) . nodeFromLine) .
        tail . B.lines <$> B.readFile nodeFl
    runResourceT $ fromLabeledEdges' edgeFl (toEdge nodeMap)
  where
    toEdge nodeMap fl = sourceFileBS fl .| linesUnboundedAsciiC .|
        (dropC 1 >> mapC f) 
      where
        f l = ( ( M.lookupDefault undefined (mk f2) nodeMap
                , M.lookupDefault undefined (mk f1) nodeMap )
              , readDouble f3 )
          where
            [f1,f2,f3,_] = B.split ',' l
{-# INLINE mkNetwork #-}

-- | Read network files as nodes and edges
readNetwork :: FilePath   -- ^ nodes
            -> FilePath   -- ^ edges
            -> IO ([NetNode], [((GeneName, GeneName), Double)])
readNetwork nodeFl edgeFl = do
    nds <- map nodeFromLine . tail . B.lines <$> B.readFile nodeFl
    es <- map f . tail . B.lines <$> B.readFile edgeFl
    return (nds, es)
  where
    f l = ( ( mk f2, mk f1), readDouble f3 )
        where
        [f1,f2,f3,_] = B.split ',' l
