-- | Get network from external source

{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network.External
    ( getExternalLinks
    ) where

import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import           Data.List.Ordered                 (nubSort)
import qualified Data.HashMap.Strict                   as M

import           Taiji.Types

-- | Get nodes and edges from the external network file.
getExternalLinks :: M.HashMap GeneName (Double, Double)  -- ^ Gene expression
                 -> FilePath                -- ^ External network file  
                 -> IO ([NetNode], [((GeneName, GeneName), Double)])
getExternalLinks expr net = do
    es <- readEdges net
    let ns = map (mkNode expr) $ nubSort $ concatMap (\(a,b) -> [a,b]) es
    return (ns, map (mkEdge expr) es)

readEdges :: FilePath -> IO [(GeneName, GeneName)]
readEdges fl = do
    c <- B.readFile fl
    return $ map (f . B.split '\t') $ B.lines c
  where
    f xs = (mk $ xs!!0, mk $ xs!!1)
{-# INLINE readEdges #-}

mkNode :: M.HashMap GeneName (Double, Double)  -- ^ Gene expression
       -> GeneName
       -> NetNode
mkNode expr nd = NetNode { _node_name = nd
                         , _node_weight = exp scaledNdExpr
                         , _node_expression = Just ndExpr }
  where
    (ndExpr, scaledNdExpr) = M.lookupDefault (0.1, 0) nd expr
{-# INLINE mkNode #-}

mkEdge :: M.HashMap GeneName (Double, Double)  -- ^ Gene expression
       -> (GeneName, GeneName)
       -> ((GeneName, GeneName), Double)
mkEdge expr (fr, to) = ((fr, to), sqrt frExpr)
  where
    (frExpr, _) = M.lookupDefault (0.1, 0) fr expr
{-# INLINE mkEdge #-}