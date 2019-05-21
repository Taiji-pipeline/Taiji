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
import qualified Taiji.SingleCell as SingleCell
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.SC.ATACSeq as SCATACSeq
import qualified Taiji.Pipeline.RNASeq as RNASeq
import qualified Taiji.Pipeline.SC.DropSeq as DropSeq

import           Taiji.Pipeline.ATACSeq.Types (ATACSeqConfig (..))
import           Taiji.Pipeline.SC.ATACSeq.Types (SCATACSeqConfig (..))
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
    _scatacseq_annotation = _taiji_annotation

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
    _dropseq_cell_barcode_length _ = 12
    _dropseq_molecular_barcode_length _ = 8

mainWith defaultMainOpts
    { programHeader = printf "Taiji-v%s" $ showVersion version } $ do
        Core.builder
        SingleCell.builder

        namespace "RNA" $ RNASeq.inputReader "RNA-seq"
        namespace "RNA" RNASeq.builder
        namespace "ATAC" ATACSeq.builder
        [ "ATAC_Get_TFBS", "ATAC_Get_Peak", "HiC_Read_Input"
            , "RNA_Make_Expr_Table" ] ~> "Create_Linkage_Prep"
        ["Create_Linkage", "RNA_Make_Expr_Table"] ~> "Compute_Ranks_Prep"

        namespace "SCATAC" SCATACSeq.builder
        namespace "DropSeq" DropSeq.builder
        [ "SCATAC_Find_TFBS", "SCATAC_Make_CutSite_Index",
            "DropSeq_Quantification" ] ~> "Compute_Ranks_SC_Prep"
        ["SCATAC_Find_TFBS", "SCATAC_Call_Peak_Cluster"] ~>
            "Compute_Ranks_SC_Cluster_Prep"