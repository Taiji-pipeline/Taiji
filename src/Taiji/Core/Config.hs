module Taiji.Core.Config where

import           Bio.Pipeline.CallPeaks
import           Bio.Pipeline.Utils
import           Control.Lens
import           Data.Default (def)
import           Taiji.Pipeline.ATACSeq.Config (ATACSeqConfig (..))
import           Taiji.Pipeline.RNASeq.Config  (RNASeqConfig (..))
import           Taiji.Types                   (TaijiConfig (..))

instance ATACSeqConfig TaijiConfig where
    _atacseq_output_dir = asDir . (++ "/ATACSeq") . _taiji_output_dir
    _atacseq_input = _taiji_input
    _atacseq_picard = _taiji_picard
    _atacseq_bwa_index = fmap (++ "/genome.fa") . _taiji_bwa_index
    _atacseq_genome_fasta = _taiji_genome
    _atacseq_genome_index = _taiji_genome_index
    _atacseq_motif_file = _taiji_motif_file
    _atacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200

instance RNASeqConfig TaijiConfig where
    _rnaseq_genome_fasta = _taiji_genome
    _rnaseq_star_index = _taiji_star_index
    _rnaseq_annotation = _taiji_annotation
    _rnaseq_rsem_index = _taiji_rsem_index
    _rnaseq_input = _taiji_input
    _rnaseq_output_dir = asDir . (++ "/RNASeq") . _taiji_output_dir
