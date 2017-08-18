{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Main where

import Control.Lens ((.=))
import           Scientific.Workflow
import           Taiji.Core.Config     (TaijiConfig)
import qualified Taiji.Pipeline.RNASeq as RNASeq
import qualified Taiji.Pipeline.ATACSeq as ATACSeq

initialization :: () -> WorkflowConfig TaijiConfig ()
initialization _ = return ()

mainWith defaultMainOpts { programHeader = "Taiji" } $ do
    namespace "RNA" RNASeq.builder
    namespace "ATAC" ATACSeq.builder
    nodeS "Initialization" 'initialization $ submitToRemote .= Just False
    ["Initialization"] ~> "RNA_Make_Index"
    ["Initialization"] ~> "ATAC_Make_Index"
