{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Extra.Functions (getTFModule) where

import           AI.Clustering.Hierarchical
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils          (asDir, getPath)
import           Control.Lens                ((^.))
import           Control.Monad
import           Control.Monad.Reader        (asks, liftIO)
import           Data.Binary                 (decodeFile)
import qualified Data.ByteString.Char8       as B
import           Data.CaseInsensitive        (CI)
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
import           Taiji.Core.Functions        (buildNet, readExpression,
                                              transform_peak_height)
import           Taiji.Types

type GeneName = CI B.ByteString

getTFModule :: ( Maybe (File '[] 'Tsv)
               , (T.Text, File '[] 'Other) )
            -> WorkflowConfig TaijiConfig (T.Text, String)
getTFModule (expr, (grp, fl)) = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        rnaseqData <- case expr of
            Nothing -> return Nothing
            Just e  -> Just <$> readExpression (e^.location) 1
        links <- decodeFile $ fl^.location
        let gr = buildNet False (B.pack $ T.unpack grp) links rnaseqData
            output = dir ++ "/" ++ T.unpack grp ++ "_network.tsv"
        return (grp, tfModule $ tfProfile gr)

tfModule :: [(GeneName, U.Vector Double)] -> String
tfModule xs = drawDendrogram $ fmap (show . fst) $ hclust Ward (V.fromList xs) (euclidean `on` snd)

tfProfile :: LGraph D (GeneName, Double) (Double, Double, Double)
          -> [(GeneName, U.Vector Double)]
tfProfile gr = mapMaybe fn $ nodes gr
  where
    fn nd | null parent = Nothing
          | otherwise = Just (fst $ nodeLab gr nd, val)
      where
        parent = pre gr nd
        val = U.create $ do
            vec <- UM.replicate n 0
            forM_ parent $ \x -> do
                let v = (\(_,a,_) -> transform_peak_height a) $ edgeLab gr (x, nd)
                UM.unsafeWrite vec x v
            return vec
    n = nNodes gr
