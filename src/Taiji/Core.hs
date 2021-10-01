{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment.Parser (mkInputReader)
import Bio.Data.Experiment.Types
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Map.Strict              as M
import Control.Workflow

import           Taiji.Core.Network
import           Taiji.Core.Ranking
import           Taiji.Core.RegulatoryElement
import Taiji.Prelude

aggregate :: MonadIO m
          => ( [ATACSeq S f1]
             , [ATACSeq S f2]
             , [HiC S f3]
             , Maybe (File '[] 'Tsv)
             , Maybe (File '[] 'Tsv) )
             -- ^ (TFBS, peaks for promoter, HiC loops, expression, estimated expression)
          -> m [ (ATACSeq S f1, f2, Maybe f3, Maybe (File '[] 'Tsv)) ]
aggregate (tfbs, atac_peaks, hic, rnaE, atacE) = liftIO $ do
    grps <- case expr of
        Nothing -> return Nothing
        Just fl -> Just . map (T.pack . B.unpack) . tail . B.split '\t' .
            head . B.lines <$> B.readFile (fl^.location)
    return $ flip mapMaybe tfbs $ \e ->
        let grp = e^.groupName._Just
            pro = M.findWithDefault undefined grp atacFileMap
            hic' = M.lookup grp hicFileMap
        in case grps of
            Nothing -> Just (e, pro, hic', expr)
            Just grps' -> if grp `elem` grps'
                then Just (e, pro, hic', expr)
                else Nothing
  where
    expr = maybe atacE Just rnaE
    atacFileMap = M.fromList $ map getFile atac_peaks
    hicFileMap = M.fromList $ map getFile hic
    getFile x = (x^.groupName._Just, x^.replicates._2.files)

builder :: Builder ()
builder = do
    node "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- mkInputReader input "HiC" (\_ x -> HiC x)
            return $ getHiCLoops hic
        |] $ doc .= "Read HiC loops from input file."

    node "Create_Linkage_Prep" 'aggregate $ return ()
    nodePar "Create_Linkage" 'saveAssociations $ memory .= 20
    path ["Create_Linkage_Prep", "Create_Linkage"]

    uNode "Compute_Ranks_Prep" [| \(x, expr) -> return $ zip x $ repeat expr |]
    nodePar "Compute_Ranks" 'computeRanks $ do
        doc .= "Perform personalized Pagerank."
        memory .= 20
    node "Output_Ranks" [| \input -> do
        dir <- asks _taiji_output_dir >>= getPath
        let output1 = dir ++ "/GeneRanks.tsv"
            output2 = dir ++ "/GeneRanks_PValues.tsv"
            output3 = dir ++ "/GeneRanks.html"
        liftIO $ outputRanks output1 output2 output3 input
        |] $ return ()
    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Ranks"]