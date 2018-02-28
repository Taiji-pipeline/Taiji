{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Extra (builder) where

import           Control.Lens
import           Scientific.Workflow

import           Taiji.Extra.Functions

builder :: Builder ()
builder = do
    nodePS 1 "TF_Module" 'getTFModule $ do
        note .= ""
    path ["TFRank_Prep", "TF_Module"]
