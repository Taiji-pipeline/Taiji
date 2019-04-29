{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser   (readHiC, readHiCTSV)
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import           Control.Lens
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Reader         (asks)
import qualified Data.Map.Strict              as M
import           Data.Maybe                   (mapMaybe)
import           Scientific.Workflow

import           Taiji.Core.Network
import           Taiji.Core.Network.SingleCell
import           Taiji.Core.Ranking
import           Taiji.Core.RegulatoryElement
import           Taiji.Types                  (_taiji_input)

aggregate :: ( [ATACSeq S f1]   -- ^ TFBS
             , [ATACSeq S f2]   -- ^ peaks for promoter
             , [HiC S f3]  -- ^ HiC loops
             , Maybe (File '[] 'Tsv) )         -- ^ Expression
          -> IO [ (ATACSeq S f1, f2, Maybe f3, Maybe (File '[] 'Tsv)) ]
aggregate (tfbs, atac_peaks, hic, expr) = do
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
    atacFileMap = M.fromList $ map getFile atac_peaks
    hicFileMap = M.fromList $ map getFile hic
    getFile x = (x^.groupName._Just, x^.replicates._2.files)

builder :: Builder ()
builder = do
    nodeS "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- if ".tsv" == reverse (take 4 $ reverse input)
                then readHiCTSV input "HiC"
                else readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ do
            submitToRemote .= Just False
            note .= "Read HiC loops from input file."

    node "Create_Linkage_Prep" 'aggregate $ do
        note .= "Prepare for parallel execution."
        submitToRemote .= Just False
    nodePS 1 "Create_Linkage" 'createLinkage $ do
        remoteParam .= "--mem=20000 -p gpu"
    path ["Create_Linkage_Prep", "Create_Linkage"]

    node' "Compute_Ranks_Prep" [| \(x, expr) -> zip x $ repeat expr |] $ do
        submitToRemote .= Just False
    nodePS 1 "Compute_Ranks" 'computeRanks $ do
        note .= "Perform personalized Pagerank."
        remoteParam .= "--mem=20000 -p gpu"
    nodeS "Output_Ranks" 'outputRanks $ return ()

    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Ranks"]


    -- Single cell
    node "Compute_Ranks_SC_Prep" 'prepDataSet $ submitToRemote .= Just False
    nodePS 1 "Compute_Ranks_SC" 'computeRanksSC $ return ()
    path ["Compute_Ranks_SC_Prep", "Compute_Ranks_SC"]

