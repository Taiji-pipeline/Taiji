{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network where

import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (asDir, getPath)
import           Bio.Utils.Functions               (scale)
import           Bio.Utils.Misc                    (readDouble)
import           Conduit
import           Control.Lens                      hiding (pre, to)
import           Control.Monad
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (CI, mk, original)
import           Data.Conduit.Cereal               (conduitGet2)
import           Data.Double.Conversion.ByteString (toShortest)
import           Data.Either                       (fromRight)
import           Data.List                         (transpose)
import           Data.List.Ordered                 (nubSort)
import qualified Data.Map.Strict                   as M
import qualified Data.Matrix.Unboxed               as MU
import           Data.Maybe                        (fromMaybe, mapMaybe)
import           Data.Monoid                       ((<>))
import           Data.Serialize                    (decode, encode, get)
import qualified Data.Text                         as T
import qualified Data.Vector.Unboxed               as U
import           IGraph
import           IGraph.Structure                  (personalizedPagerank)
import           Scientific.Workflow               hiding (_data)

import           Taiji.Core.Config                 ()
import           Taiji.Types

computeRanks :: ( ATACSeq S (File '[] 'Other)
                , Maybe (File '[] 'Tsv) )          -- ^ Expression
             -> WorkflowConfig TaijiConfig (T.Text, File '[] 'Other)
computeRanks (atac, expr) = do
    dir <- asks (asDir . _taiji_output_dir) >>= getPath . (<> asDir "/Network")
    liftIO $ do
        network <- mkNetwork $ runIdentity (atac^.replicates) ^. files.location
        gr <- fmap pageRank $ case expr of
            Nothing -> return network
            Just e -> do
                expr' <- readExpression 1
                    (B.pack $ T.unpack grp ) $ e^.location
                return $ assignWeights expr' network
        let output = dir ++ "/" ++ T.unpack grp ++
                "_network.bin"
        B.writeFile output $ encode gr
        return (grp, location .~ output $ emptyFile)
  where
    grp = atac^.groupName._Just

mkNetwork :: FilePath -> IO (Graph 'D NetNode NetEdge)
mkNetwork input = runResourceT $ fromLabeledEdges' input toEdge
  where
    toEdge fl = sourceFileBS fl .| conduitGet2 get .| concatMapC f
      where
        f (gene, tfs1, tfs2) = zipWith g (repeat gene) $ M.toList $ M.fromListWith max $ tfs1 ++ tfs2
        g a b = ( (NetNode a Nothing Nothing Nothing, NetNode (fst b) Nothing Nothing Nothing)
                , NetEdge Nothing Nothing (snd b) )
{-# INLINE mkNetwork #-}

assignWeights :: M.Map (CI B.ByteString) (Double, Double)   -- ^ Gene expression
              -> Graph 'D NetNode NetEdge
              -> Graph 'D NetNode NetEdge
assignWeights weights gr = emap assignEdgeWeight $
    nmap assignNodeWeight gr
  where
    assignNodeWeight (_, x) =
        let (raw, scaled) = M.findWithDefault (0.1, -10) (nodeName x) weights
        in x{nodeExpression=Just raw, nodeScaledExpression=Just scaled}
    assignEdgeWeight ((_, to), x) =
        let e = M.findWithDefault (0.1, -10) (nodeName $ nodeLab gr to) weights
        in x{weightExpression=Just $ fst e}
{-# INLINE assignWeights #-}

-- | Read RNA expression data
readExpression :: Double    -- ^ Threshold to call a gene as non-expressed
               -> B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.Map (CI B.ByteString) (Double, Double)) -- ^ absolute value and z-score
readExpression cutoff ct fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (U.fromList . map ((+pseudoCount) . readDouble) . tail) dat
        dataTable' = MU.fromRows $ zipWith U.zip dataTable $
            map computeZscore dataTable :: MU.Matrix (Double, Double)
        idx = fromMaybe (error $ "Cell type:" ++ B.unpack ct ++ " not found!") $
            lookup ct $ zip header [0..]
    return $ M.fromList $ zip rowNames $ U.toList $ dataTable' `MU.takeColumn` idx
  where
    pseudoCount = 0.1
    computeZscore xs
        | U.all (<cutoff) xs || U.all (== U.head xs) xs = U.replicate (U.length xs) (-10)
        | otherwise = scale xs
{-# INLINE readExpression #-}

pageRank :: Graph 'D NetNode NetEdge
         -> Graph 'D NetNode NetEdge
pageRank gr = nmap (\(i, x) -> x{pageRankScore=Just $ ranks U.! i}) gr
  where
    labs = map (nodeLab gr) $ nodes gr
    nodeWeights = map (transform_node_weight . fromMaybe (-10) .
        nodeScaledExpression) labs
    edgeWeights = map (combine . edgeLab gr) $ edges gr
    ranks = U.fromList $ personalizedPagerank gr nodeWeights (Just edgeWeights) 0.85
    combine NetEdge{..} = transform_exp (fromMaybe 1 weightExpression) *
        getSiteWeight sites
{-# INLINE pageRank #-}


outputRanks :: [(T.Text, File '[] 'Other)]
            -> WorkflowConfig TaijiConfig FilePath
outputRanks inputs = do
    dir <- asks _taiji_output_dir >>= getPath . asDir
    let output = dir ++ "/GeneRanks.tsv"

    ranks <- forM inputs $ \(_, fl) -> liftIO $ do
        gr <- fmap (fromRight (error "decode fail") . decode) $
            B.readFile $ fl^.location :: IO (Graph 'D NetNode NetEdge)
        return $! M.fromList $ flip mapMaybe (nodes gr) $ \i ->
            if null (pre gr i)
                then Nothing
                else do
                    let n = nodeLab gr i
                    x <- pageRankScore n
                    return (nodeName n, x)
    let genes = nubSort $ concatMap M.keys ranks
        header = B.pack $ T.unpack $ T.intercalate "\t" $
            "Gene" : fst (unzip inputs)
        ranks' = flip map ranks $ \r -> flip map genes $
            \g -> M.findWithDefault 0 g r

    liftIO $ B.writeFile output $ B.unlines $
        header : zipWith toBS (map original genes) (transpose ranks')
    return output
  where
    toBS nm xs = B.intercalate "\t" $ nm : map toShortest xs

--------------------------------------------------------------------------------
-- Transformation functions
--------------------------------------------------------------------------------

getSiteWeight :: Double -> Double
getSiteWeight = id
{-# INLINE getSiteWeight #-}

transform_exp :: Double -> Double
transform_exp = sqrt
{-# INLINE transform_exp #-}

transform_corr :: Double -> Double
transform_corr = abs
{-# INLINE transform_corr #-}

transform_node_weight :: Double -> Double
transform_node_weight = exp
{-# INLINE transform_node_weight #-}
