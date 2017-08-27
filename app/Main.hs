{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Main where

import           Control.Lens           ((.=))
import           Scientific.Workflow
import           Taiji.Core             (builder)
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.RNASeq  as RNASeq
import           Taiji.Types            (TaijiConfig)

initialization :: () -> WorkflowConfig TaijiConfig ()
initialization _ = return ()

mainWith defaultMainOpts { programHeader = "Taiji" } $ do
    namespace "RNA" RNASeq.builder
    namespace "ATAC" ATACSeq.builder
    builder
    nodeS "Initialization" 'initialization $ submitToRemote .= Just False
    ["Initialization"] ~> "RNA_Read_Input"
    ["Initialization"] ~> "ATAC_Read_Input"
    path ["ATAC_Merge_Bed", "Find_Active_Promoter"]
    ["Find_Active_Promoter", "ATAC_Get_TFBS"] ~> "Link_Gene_TF_Prep"
    ["RNA_Make_Expr_Table", "Link_Gene_TF"] ~> "TFRank_Prep"
    ["Output_Rank", "RNA_Make_Expr_Table", "Link_Gene_TF"] ~> "Export_Results"
