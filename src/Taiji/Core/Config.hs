{-# LANGUAGE DeriveGeneric #-}
module Taiji.Core.Config where

import           Bio.Pipeline.Utils
import           Data.Aeson
import Data.Aeson.Types (fieldLabelModifier)
import           Data.Default                  (Default (..))
import           GHC.Generics                  (Generic)
import           Taiji.Pipeline.ATACSeq.Config (ATACSeqConfig (..))
import           Taiji.Pipeline.RNASeq.Config  (RNASeqConfig (..))

data TaijiConfig = TaijiConfig
    { _taiji_output_dir :: Directory
    , _taiji_input      :: FilePath
    , _taiji_picard     :: Maybe FilePath
    , _taiji_genome     :: Maybe FilePath
    , _taiji_bwa_index  :: Maybe FilePath
    , _taiji_star_index :: Maybe FilePath
    , _taiji_annotation :: Maybe FilePath
    , _taiji_rsem_index :: Maybe FilePath
    } deriving (Generic)

instance Default TaijiConfig where
    def = TaijiConfig
        { _taiji_output_dir = asDir "output"
        , _taiji_input = "input.yml"
        , _taiji_picard = def
        , _taiji_genome = def
        , _taiji_bwa_index = def
        , _taiji_star_index = def
        , _taiji_annotation = def
        , _taiji_rsem_index = def
        }

instance ToJSON TaijiConfig where
    toEncoding = genericToEncoding defaultOptions
        { fieldLabelModifier = drop 7 }

instance FromJSON TaijiConfig where
    parseJSON = genericParseJSON defaultOptions
        { fieldLabelModifier = ("_taiji_" ++) }

instance ATACSeqConfig TaijiConfig where
    _atacseq_output_dir = _taiji_output_dir
    _atacseq_input = _taiji_input
    _atacseq_picard = _taiji_picard
    _atacseq_bwa_index = _taiji_bwa_index
    _atacseq_genome_fasta = _taiji_genome

instance RNASeqConfig TaijiConfig where
    _rnaseq_genome_fasta = _taiji_genome
    _rnaseq_star_index = _taiji_star_index
    _rnaseq_annotation = _taiji_annotation
    _rnaseq_rsem_index = _taiji_rsem_index
    _rnaseq_input = _taiji_input
    _rnaseq_output_dir = _taiji_output_dir
