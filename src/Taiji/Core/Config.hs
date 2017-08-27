{-# LANGUAGE DeriveGeneric #-}
module Taiji.Core.Config where

import           Bio.Pipeline.Utils
import           Data.Aeson
import Data.Aeson.Types (fieldLabelModifier)
import           Data.Default                  (Default (..))
import           GHC.Generics                  (Generic)
import           Taiji.Pipeline.ATACSeq.Config (ATACSeqConfig (..))
import           Taiji.Pipeline.RNASeq.Config  (RNASeqConfig (..))
import Taiji.Types (TaijiConfig(..))

instance ATACSeqConfig TaijiConfig where
    _atacseq_output_dir = asDir . (++ "/ATACSeq") . _taiji_output_dir
    _atacseq_input = _taiji_input
    _atacseq_picard = _taiji_picard
    _atacseq_bwa_index = _taiji_bwa_index
    _atacseq_genome_fasta = _taiji_genome
    _atacseq_genome_index = _taiji_genome_index
    _atacseq_motif_file = _taiji_motif_file

instance RNASeqConfig TaijiConfig where
    _rnaseq_genome_fasta = _taiji_genome
    _rnaseq_star_index = _taiji_star_index
    _rnaseq_annotation = _taiji_annotation
    _rnaseq_rsem_index = _taiji_rsem_index
    _rnaseq_input = _taiji_input
    _rnaseq_output_dir = asDir . (++ "/RNASeq") . _taiji_output_dir
