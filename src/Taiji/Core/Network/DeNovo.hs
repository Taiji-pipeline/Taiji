-- | Infer network de novo from data

{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network.DeNovo
    ( createLinkage
    ) where

import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (getPath, asDir)
import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Conduit
import qualified Data.IntervalMap.Strict as IM
import           Control.Lens                      hiding (pre, to)
import qualified Data.Set as S
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import Data.Ord (comparing)
import           Data.List                         (foldl', maximumBy)
import qualified Data.Map.Strict                   as M
import           Data.Maybe                        (fromJust, mapMaybe)
import qualified Data.Text                         as T
import           Scientific.Workflow               hiding (_data)
import System.IO

import Taiji.Core.Network.Utils
import Taiji.Core.RegulatoryElement (findTargets)
import           Taiji.Core.Config                 ()
import           Taiji.Types

createLinkage :: ( ATACSeq S ( File tag1 'Bed         -- ^ Active promoters
                             , File tag2 'Bed         -- ^ TFBS
                             )
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ promoter activity
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ enhancer activity
                 , Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
                 , Maybe (File '[] 'Tsv)          -- ^ Expression
                 )
              -> WorkflowConfig TaijiConfig
                    (ATACSeq S (File '[] 'Other, File '[] 'Other))
createLinkage (atac, pro, enh, hic, expr) = do
    dir <- asks
        ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir)
        >>= getPath
    let netEdges = dir ++ "/edges_combined.csv"
        netNodes = dir ++ "/nodes.csv"
        bindingEdges = dir ++ "/edges_binding.csv"
    liftIO $ do
        expr' <- case expr of
            Nothing -> return M.empty
            Just e -> readExpression 1 (B.pack $ T.unpack grp ) $ e^.location
        activityPro <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_pro
        activityEnh <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_enh
        withFile netEdges WriteMode $ \h1 -> withFile bindingEdges WriteMode $ \h2 -> do
            B.hPutStrLn h1 ":START_ID,:END_ID,weight,:TYPE"
            B.hPutStrLn h2 $ ":START_ID,:END_ID,chr,start:int,end:int," <>
                "annotation,affinity,:TYPE"
            let conduit = findTargets active_pro tfbs hic .|
                    createLinkage_ activityPro activityEnh expr' .|
                    mapM_C (liftIO . outputEdge h1 h2)
            s <- execStateT (runConduit conduit) S.empty
            let nodeHeader = "geneName:ID,expression,expressionZScore"
            B.writeFile netNodes $ B.unlines $ (nodeHeader:) $
                map nodeToLine $ S.toList s
    return $ atac & replicates.mapped.files .~
        ( emptyFile & location .~ netNodes
        , emptyFile & location .~ netEdges )
  where
    outputEdge h1 h2 e = B.hPutStrLn hdl $ edgeToLine e
      where
        hdl = case _edge_type e of
            Combined _ -> h1
            Binding{..} -> h2
    fl_pro = either (^.location) (^.location) pro
    fl_enh = either (^.location) (^.location) enh
    (active_pro, tfbs) = atac^.replicates._2.files
    grp = atac^.groupName._Just
{-# INLINE createLinkage #-}

createLinkage_ :: BEDTree Double   -- ^ Promoter activities
               -> BEDTree Double   -- ^ Enhancer activities
               -> M.Map GeneName (Double, Double)   -- ^ Gene expression
               -> ConduitT (GeneName, ([BED], [BED]))
                           NetEdge
                           (StateT (S.Set NetNode) IO) ()
createLinkage_ act_pro act_enh expr = concatMapMC $ \(geneName, (ps, es)) -> do
    let tfEnhancer = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (g act_enh) es
        edgeEnhancer = flip concatMap tfEnhancer $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "enhancer"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfPromoter = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (g act_pro) ps
        edgePromoter = flip concatMap tfPromoter $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "promoter"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfs = M.toList $ fmap (lp 2 . map (fromJust . (^.score))) $
            M.fromListWith (++) $ tfEnhancer ++ tfPromoter
        (geneExpr, scaledGeneExpr) = M.findWithDefault (0.1, 0) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_expression = Just geneExpr
                           , _node_scaled_expression = Just scaledGeneExpr }
    modify' $ S.insert geneNode
    edgeCombined <- forM tfs $ \(tfName, w) -> do
        let (tfExpr, scaledTfExpr) = M.findWithDefault (0.1, 0) tfName expr
            tfNode = NetNode { _node_name = tfName
                             , _node_expression = Just tfExpr
                             , _node_scaled_expression = Just scaledTfExpr }
        modify' $ S.insert tfNode
        return $ NetEdge { _edge_from = tfName
                         , _edge_to = geneName
                         , _edge_type = Combined (w * sqrt tfExpr) }
    return $ edgePromoter ++ edgeEnhancer ++ edgeCombined
  where
    getBestMotif xs = runIdentity $ runConduit $
        mergeBedWith (maximumBy (comparing (^.score))) xs .| sinkList
    g act bed = case IM.elems (intersecting act bed) of
        [] -> Nothing
        xs -> Just ( getTFName bed
            , [bed & score.mapped %~ (transform_peak_height (maximum xs) *) . transform_site_pvalue] )
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE createLinkage_ #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}

