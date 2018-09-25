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
import qualified Data.Map.Strict                   as M

import           Taiji.Core.Config                 ()
import           Taiji.Types

-- | Get nodes and edges from the external network file.
getExternalLinks :: M.Map GeneName (Double, Double)  -- ^ Gene expression
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
    f [a,b] = (mk a, mk b)
    f _ = undefined
{-# INLINE readEdges #-}

mkNode :: M.Map GeneName (Double, Double)  -- ^ Gene expression
       -> GeneName
       -> NetNode
mkNode expr nd = NetNode { _node_name = nd
                         , _node_expression = Just ndExpr
                         , _node_scaled_expression = Just scaledNdExpr }
  where
    (ndExpr, scaledNdExpr) = M.findWithDefault (0.1, 0) nd expr
{-# INLINE mkNode #-}

mkEdge :: M.Map GeneName (Double, Double)  -- ^ Gene expression
       -> (GeneName, GeneName)
       -> ((GeneName, GeneName), Double)
mkEdge expr (fr, to) = ((fr, to), frExpr)
  where
    (frExpr, _) = M.findWithDefault (0.1, 0) fr expr
{-# INLINE mkEdge #-}