{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Types where

import           Bio.Data.Bed
import           Bio.Pipeline.Instances ()
import           Control.Lens           (makeLenses)
import           Data.Aeson
import IGraph
import qualified Data.ByteString.Char8  as B
import           Data.CaseInsensitive   (CI)
import           Data.Default.Class
import           Data.Hashable
import qualified Data.Map.Strict        as M
import qualified Data.Matrix.Unboxed    as MU
import           Data.Serialize         (Serialize (..))
import           Data.Serialize.Text    ()
import qualified Data.Text              as T
import qualified Data.Vector            as V
import           Data.Vector.Serialize  ()
import           GHC.Generics           (Generic)

type Promoter = BEDExt BED3 (CI B.ByteString)
type RegDomain = BEDExt BED3 (CI B.ByteString)

makeLenses ''SiteInfo

data RankTable = RankTable
    { rowNames    :: V.Vector T.Text
    , colNames    :: V.Vector T.Text
    , ranks       :: MU.Matrix Double
    , expressions :: Maybe (MU.Matrix Double)
    } deriving (Generic)

data TaijiResults = TaijiResults
    { ranktable :: RankTable
    , nets      :: M.Map T.Text (M.Map T.Text T.Text)
    } deriving (Generic)

data TaijiConfig = TaijiConfig
    { _taiji_output_dir   :: FilePath
    , _taiji_input        :: FilePath
    , _taiji_picard       :: Maybe FilePath
    , _taiji_genome       :: Maybe FilePath
    , _taiji_bwa_index    :: Maybe FilePath
    , _taiji_star_index   :: Maybe FilePath
    , _taiji_annotation   :: Maybe FilePath
    , _taiji_rsem_index   :: Maybe FilePath
    , _taiji_genome_index :: Maybe FilePath
    , _taiji_motif_file   :: Maybe FilePath
    } deriving (Generic)

instance Default TaijiConfig where
    def = TaijiConfig
        { _taiji_output_dir = "output"
        , _taiji_input = "input.yml"
        , _taiji_picard = def
        , _taiji_genome = def
        , _taiji_bwa_index = def
        , _taiji_star_index = def
        , _taiji_annotation = def
        , _taiji_rsem_index = def
        , _taiji_genome_index = def
        , _taiji_motif_file = def
        }

instance ToJSON TaijiConfig where
    toJSON = genericToJSON defaultOptions
        { fieldLabelModifier = drop 7 }
    toEncoding = genericToEncoding defaultOptions
        { fieldLabelModifier = drop 7 }

instance FromJSON TaijiConfig where
    parseJSON = genericParseJSON defaultOptions
        { fieldLabelModifier = drop 7 }

data NetNode = NetNode
    { nodeName             :: CI B.ByteString
    , nodeExpression       :: Maybe Double
    , nodeScaledExpression :: Maybe Double
    , pageRankScore        :: Maybe Double
    } deriving (Generic, Show, Read, Eq)

instance Hashable NetNode where
    hashWithSalt salt at = hashWithSalt salt $ nodeName at

data NetEdge = NetEdge
    { weightExpression  :: Maybe Double
    , weightCorrelation :: Maybe Double
    , sites             :: (Maybe (Double, Int), Maybe (Double, Int), Maybe (Double, Int))
    } deriving (Generic, Show, Read)

instance Serialize NetNode
instance Serialize NetEdge
instance Serialize TaijiResults
instance Serialize RankTable

instance FromJSON TaijiResults
instance ToJSON TaijiResults
instance FromJSON RankTable
instance ToJSON RankTable

instance Default (CI B.ByteString) where
    def = ""

instance Graph d => Serialize (LGraph d NetNode NetEdge) where
    put gr = do
        put nlabs
        put es
        put elabs
      where
        nlabs = map (nodeLab gr) $ nodes gr
        es = edges gr
        elabs = map (edgeLab gr) es
    get = do
        nlabs <- get
        es <- get
        elabs <- get
        return $ mkGraph nlabs $ zip es elabs
