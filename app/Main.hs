{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Main where

import           Control.Lens           ((.=))
import           Data.Version           (showVersion)
import           Paths_Taiji            (version)
import           Scientific.Workflow
import           Text.Printf            (printf)

import           Taiji.Core             (builder)
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.RNASeq  as RNASeq
import           Taiji.Types            (TaijiConfig)

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" $ showVersion version } $ do
        namespace "RNA" RNASeq.builder
        namespace "ATAC" ATACSeq.builder
        builder
        path ["ATAC_Call_Peak", "Find_Active_Promoter"]
        ["Find_Active_Promoter", "ATAC_Get_TFBS"] ~> "Link_Gene_TF_Prep"
        ["RNA_Make_Expr_Table", "Link_Gene_TF"] ~> "TFRank_Prep"
        ["Output_Rank", "RNA_Make_Expr_Table", "Link_Gene_TF"] ~> "Export_Results"
