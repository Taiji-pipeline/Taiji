{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where

import           Control.Lens           ((.=))
import           Data.Version           (showVersion)
import           Paths_Taiji            (version)
import           Scientific.Workflow
import           Text.Printf            (printf)

import qualified Taiji.Core             as Core
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.ChIPSeq as ChIPSeq
import qualified Taiji.Pipeline.RNASeq.Classic  as RNASeq
import           Taiji.Types            (TaijiConfig)

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" $ showVersion version } $ do
        namespace "RNA" $ RNASeq.inputReader "RNA-seq"
        namespace "RNA" RNASeq.builder
        namespace "ATAC" ATACSeq.builder
        namespace "H3K27ac" $ ChIPSeq.inputReader "H3K27ac"
        namespace "H3K27ac" ChIPSeq.builder
        Core.builder
        path ["ATAC_Call_Peak", "Find_Active_Promoter"]
        [ "Find_Active_Promoter", "ATAC_Get_TFBS", "ATAC_Call_Peak"
            , "H3K27ac_Get_Peak", "HiC_Read_Input", "RNA_Make_Expr_Table"
            ] ~> "Create_Linkage_Prep"
        ["Create_Linkage", "RNA_Make_Expr_Table"] ~> "Compute_Ranks_Prep"
