{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Extra (builder) where

import           Control.Lens
import           Scientific.Workflow

import           Taiji.Extra.Functions

builder :: Builder ()
builder = do
    nodePS 1 "TF_Module" 'getTFModule $ note .= ""
    path ["Compute_Ranks", "TF_Module"]
