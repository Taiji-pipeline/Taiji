{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
module Taiji.Extra.Functions
    ( getTFModule
    ) where

import           AI.Clustering.Hierarchical
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (asDir, getPath)
import           Control.Lens                      ((^.))
import           Control.Monad
import           Bio.Utils.Misc                    (readDouble)
import           Control.Monad.Reader              (asks, liftIO)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI)
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Either                       (fromRight)
import           Data.Function                     (on)
import qualified Data.HashMap.Strict               as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe
import           Data.Monoid                       ((<>))
import           Data.Serialize                    (decode)
import qualified Data.Text                         as T
import qualified Data.Vector                       as V
import qualified Data.Vector.Unboxed               as U
import qualified Data.Vector.Unboxed.Mutable       as UM
import           IGraph
import           Scientific.Workflow

import           Taiji.Core.Config                 ()
import           Taiji.Core.Network                (getSiteWeight)
import           Taiji.Types

type GeneName = CI B.ByteString

getTFModule :: (T.Text, File '[] 'Other)
            -> WorkflowConfig TaijiConfig (T.Text, String)
getTFModule (grp, fl) = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        gr <- fmap (fromRight undefined . decode) $ B.readFile $
            fl^.location :: IO (LGraph D NetNode NetEdge)
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

-- | Read RNA expression data
readExpression :: B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.HashMap GeneName Double)
readExpression ct fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = MU.fromRows $ map (U.fromList . map readDouble . tail) dat
            :: MU.Matrix Double
        idx = fromMaybe (error $ "Cell type:" ++ B.unpack ct ++ " not found!") $
            lookup ct $ zip header [0..]
    return $ M.fromList $ zip rowNames $ U.toList $ dataTable `MU.takeColumn` idx
{-# INLINE readExpression #-}
