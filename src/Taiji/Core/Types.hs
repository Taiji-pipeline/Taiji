{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE FlexibleInstances #-}

module Taiji.Core.Types where

import           Bio.Data.Bed          (BED)
import           Data.Binary           (Binary (..))
import qualified Data.ByteString.Char8 as B
import           Data.CaseInsensitive  (CI, mk, original)
import qualified Data.Map.Strict       as M
import qualified Data.Matrix.Unboxed   as MU
import qualified Data.Text             as T
import           Data.Vector.Binary    ()
import qualified Data.Vector as V
import           GHC.Generics          (Generic)

type GeneName = CI B.ByteString

instance Binary (CI B.ByteString) where
    put = put . original
    get = fmap mk get

-- | Gene and its regulators
type Linkage = (GeneName, [(GeneName, [BED])])

data RankTable = RankTable
    { rowNames    :: V.Vector B.ByteString
    , colNames    :: V.Vector T.Text
    , ranks       :: MU.Matrix Double
    , expressions :: MU.Matrix Double
    } deriving (Generic)

instance Binary (MU.Matrix Double)
instance Binary RankTable

data TaijiResults = TaijiResults
    { ranktable :: RankTable
    , nets      :: M.Map T.Text (M.Map B.ByteString B.ByteString)
    } deriving (Generic)

instance Binary TaijiResults
