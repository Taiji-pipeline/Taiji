module Taiji.Prelude
    ( module Bio.Data.Experiment
    , module Conduit
    , module Data.Maybe
    , module Data.List
    , module Scientific.Workflow
    , module Taiji.Types
    , edge_weight_cutoff
    , (^.)
    , (.~)
    , (&)
    , (%%~)
    , (%~)
    , (.=)
    , mapped
    , traverse
    , _1
    , _2
    , _3
    , _Just
    ) where

import           Bio.Data.Experiment
import           Conduit
import           Control.Lens
import Data.List
import           Data.Maybe
import           Scientific.Workflow
import           Taiji.Types

-- | Cutoff for edge weights
edge_weight_cutoff :: Double
edge_weight_cutoff = 0.2