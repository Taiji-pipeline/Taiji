{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where

import           Data.Version           (showVersion)
import           Paths_Taiji            (version)
import           Scientific.Workflow
import           Text.Printf            (printf)
import           Bio.Pipeline.CallPeaks
import           Control.Lens
import           Data.Default                         (def)
import Data.Maybe

import qualified Taiji.Core             as Core
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.SC.ATACSeq as SCATACSeq
import qualified Taiji.Pipeline.ChIPSeq as ChIPSeq
import qualified Taiji.Pipeline.RNASeq as RNASeq
import qualified Taiji.Pipeline.SC.DropSeq as DropSeq

import           Taiji.Pipeline.ATACSeq.Types (ATACSeqConfig (..))
import           Taiji.Pipeline.SC.ATACSeq.Types (SCATACSeqConfig (..))
import           Taiji.Pipeline.ChIPSeq.Config        (ChIPSeqConfig (..))
import           Taiji.Pipeline.RNASeq.Config (RNASeqConfig (..))
import           Taiji.Pipeline.SC.DropSeq.Types (DropSeqConfig (..))

import           Taiji.Types                          (TaijiConfig (..))

instance ATACSeqConfig TaijiConfig where
    _atacseq_output_dir = (<> "/ATACSeq") . _taiji_output_dir
    _atacseq_input = _taiji_input
    _atacseq_bwa_index = fmap (++ "/genome.fa") . _taiji_bwa_index
    _atacseq_genome_fasta = _taiji_genome
    _atacseq_genome_index = _taiji_genome_index
    _atacseq_motif_file = _taiji_motif_file
    _atacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
                                   & cutoff .~ PValue 0.01
                                   & callSummits .~ True

instance SCATACSeqConfig TaijiConfig where
    _scatacseq_output_dir = (<> "/SCATACSeq") . _taiji_output_dir
    _scatacseq_input = _taiji_input
    _scatacseq_bwa_index = fmap (++ "/genome.fa") . _taiji_bwa_index
    _scatacseq_genome_fasta = _taiji_genome
    _scatacseq_genome_index = _taiji_genome_index
    _scatacseq_motif_file = _taiji_motif_file
    _scatacseq_callpeak_opts _ = def & mode .~ NoModel (-100) 200
                                   & cutoff .~ PValue 0.01
                                   & callSummits .~ True

instance ChIPSeqConfig TaijiConfig where
    _chipseq_output_dir = (<> "/ChIPSeq") . _taiji_output_dir
    _chipseq_input = _taiji_input
    _chipseq_bwa_index = fmap (++ "/genome.fa") . _taiji_bwa_index
    _chipseq_genome_fasta = _taiji_genome
    _chipseq_genome_index = _taiji_genome_index

instance RNASeqConfig TaijiConfig where
    _rnaseq_genome_fasta = _taiji_genome
    _rnaseq_star_index = _taiji_star_index
    _rnaseq_annotation = _taiji_annotation
    _rnaseq_rsem_index = fmap (++ "/genome") . _taiji_rsem_index
    _rnaseq_input = _taiji_input
    _rnaseq_output_dir = (<> "/RNASeq") . _taiji_output_dir

instance DropSeqConfig TaijiConfig where
    _dropseq_genome_fasta = _taiji_genome
    _dropseq_star_index = fromJust . _taiji_star_index
    _dropseq_annotation = fromJust . _taiji_annotation
    _dropseq_input = _taiji_input
    _dropseq_output_dir = (<> "/DropSeq") . _taiji_output_dir

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" $ showVersion version } $ do
        namespace "RNA" $ RNASeq.inputReader "RNA-seq"
        namespace "RNA" RNASeq.builder
        namespace "DropSeq" DropSeq.builder
        namespace "ATAC" ATACSeq.builder
        namespace "SCATAC" SCATACSeq.builder
        namespace "H3K27ac" $ ChIPSeq.inputReader "H3K27ac"
        namespace "H3K27ac" ChIPSeq.builder
        Core.builder
        path ["ATAC_Get_Peak", "Find_Active_Promoter"]
        [ "Find_Active_Promoter", "ATAC_Get_TFBS", "ATAC_Get_Peak"
            , "H3K27ac_Get_Peak", "HiC_Read_Input", "RNA_Make_Expr_Table"
            ] ~> "Create_Linkage_Prep"
        ["Create_Linkage", "RNA_Make_Expr_Table"] ~> "Compute_Ranks_Prep"
