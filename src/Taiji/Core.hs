{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}
module Taiji.Core (builder) where

import Scientific.Workflow

import Taiji.Core.Functions

builder :: Builder ()
builder = do
    nodePS 1 "Find_Active_Promoter" 'getActivePromoter $ return ()
