{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Extra.Functions
    ( getTFModule
    ) where

import           AI.Clustering.Hierarchical
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils          (asDir, getPath)
import           Conduit
import           Control.Lens                ((^.))
import           Control.Monad
import           Control.Monad.Reader        (asks, liftIO)
import           Data.Function               (on)
import           Data.Maybe
import           Data.Monoid                 ((<>))
import qualified Data.Text                   as T
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Mutable as UM
import           IGraph
import           Scientific.Workflow

import           Taiji.Core.Config           ()
import           Taiji.Core.Network          (getSiteWeight)
import           Taiji.Types

getTFModule :: (T.Text, File '[] 'Other)
            -> WorkflowConfig TaijiConfig (T.Text, String)
getTFModule (grp, fl) = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        gr <- runResourceT $ runConduit $ sourceFileBS (fl^.location) .|
            decodeC :: IO (LGraph D NetNode NetEdge)
        return (grp, tfModule $ tfProfile gr)

tfModule :: [(GeneName, U.Vector Double)] -> String
tfModule xs = drawDendrogram $ fmap (show . fst) $ hclust Ward (V.fromList xs) (euclidean `on` snd)

-- | TF's regulatees
tfProfile :: LGraph D NetNode NetEdge
          -> [(GeneName, U.Vector Double)]
tfProfile gr = mapMaybe fn $ nodes gr
  where
    fn nd | null parent = Nothing
          | otherwise = Just (nodeName $ nodeLab gr nd, val)
      where
        parent = pre gr nd
        val = U.create $ do
            vec <- UM.replicate n 0
            forM_ parent $ \x -> do
                let v = getSiteWeight $ sites $ edgeLab gr (x, nd)
                UM.unsafeWrite vec x v
            return vec
    n = nNodes gr
