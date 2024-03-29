{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where

import           Data.Version           (showVersion)
import           Bio.Pipeline.CallPeaks
import Data.Yaml (decodeFileThrow, Value(..))
import           Data.Default                         (def)
import qualified Data.HashMap.Strict as M
import Control.Workflow.Coordinator.Remote (Remote, RemoteConfig(..), getDefaultRemoteConfig)
import           Control.Workflow
import Control.Workflow.Main
import Data.Proxy (Proxy(..))
import qualified Data.Text as T

import           Paths_Taiji            (version)
import qualified Taiji.Core             as Core
import qualified Taiji.SingleCell as SingleCell
import qualified Taiji.Pipeline.ATACSeq as ATACSeq
import qualified Taiji.Pipeline.SC.ATACSeq as SCATACSeq
import qualified Taiji.Pipeline.RNASeq as RNASeq
import qualified Taiji.Pipeline.SC.RNASeq as SCRNASeq
import qualified Taiji.Pipeline.Fusion as Fusion

import           Taiji.Prelude
import           Taiji.Pipeline.ATACSeq.Types (ATACSeqConfig (..))
import           Taiji.Pipeline.SC.ATACSeq.Types (SCATACSeqConfig (..))
import           Taiji.Pipeline.RNASeq.Types (RNASeqConfig (..))
import           Taiji.Pipeline.SC.RNASeq.Types (SCRNASeqConfig (..))

instance ATACSeqConfig TaijiConfig where
    _atacseq_output_dir = (<> "/ATACSeq") . _taiji_output_dir
    _atacseq_input = _taiji_input
    _atacseq_assembly = _taiji_assembly
    _atacseq_bwa_index = Just . _taiji_bwa_index
    _atacseq_genome_fasta = _taiji_genome
    _atacseq_genome_index = Just . _taiji_genome_index
    _atacseq_motif_file = _taiji_motif_file
    _atacseq_bwa_seed_length = fromMaybe 32 . _taiji_bwa_seed_length
    _atacseq_callpeak_opts config = def & mode .~ NoModel (-100) 200
        & cutoff .~ QValue (fromMaybe 0.05 $ _taiji_callpeak_fdr config)
        & callSummits .~ True
        & gSize .~ _taiji_callpeak_genome_size config
        & tmpDir .~ fromMaybe "./" (_atacseq_tmp_dir config)
    _atacseq_annotation = _taiji_annotation
    _atacseq_tmp_dir = _taiji_tmp_dir

instance SCATACSeqConfig TaijiConfig where
    _scatacseq_output_dir = (<> "/SCATACSeq") . _taiji_output_dir
    _scatacseq_input = _taiji_input
    _scatacseq_batch_info = _taiji_batch_info
    _scatacseq_assembly = _taiji_assembly
    _scatacseq_bwa_index = Just . _taiji_bwa_index
    _scatacseq_genome_fasta = _taiji_genome
    _scatacseq_genome_index = Just . _taiji_genome_index
    _scatacseq_motif_file = _taiji_motif_file
    _scatacseq_callpeak_opts config = def & mode .~ NoModel (-100) 200
        & cutoff .~ QValue (fromMaybe 0.01 $ _taiji_callpeak_fdr config)
        & callSummits .~ True
        & gSize .~ _taiji_callpeak_genome_size config
        & tmpDir .~ fromMaybe "./" (_scatacseq_tmp_dir config)
    _scatacseq_annotation = _taiji_annotation
    _scatacseq_tmp_dir = _taiji_tmp_dir
    _scatacseq_blacklist = _taiji_blacklist
    _scatacseq_te_cutoff = _scatac_te_cutoff . _taiji_scatac_options
    _scatacseq_minimal_fragment = _scatac_minimal_fragment . _taiji_scatac_options
    _scatacseq_cluster_resolution_list = _scatac_cluster_resolution_list . _taiji_scatac_options
    _scatacseq_cluster_resolution = _scatac_cluster_resolution . _taiji_scatac_options
    _scatacseq_do_subclustering = _scatac_do_subclustering . _taiji_scatac_options
    _scatacseq_subcluster_resolution = _scatac_subcluster_resolution . _taiji_scatac_options
    _scatacseq_cluster_optimizer = _scatac_cluster_optimizer . _taiji_scatac_options
    _scatacseq_doublet_score_cutoff = _scatac_doublet_score_cutoff . _taiji_scatac_options
    _scatacseq_cluster_by_window = _scatac_cluster_by_window . _taiji_scatac_options
    _scatacseq_cluster_exclude = _scatac_cluster_exclude . _taiji_scatac_options
    _scatacseq_window_size = _scatac_window_size . _taiji_scatac_options
    _scatacseq_cell_barcode_length = _scatac_cell_barcode_length . _taiji_scatac_options

instance RNASeqConfig TaijiConfig where
    _rnaseq_assembly = _taiji_assembly
    _rnaseq_genome_fasta = _taiji_genome
    _rnaseq_star_index = Just . _taiji_star_index
    _rnaseq_annotation = _taiji_annotation
    _rnaseq_rsem_index = Just . _taiji_rsem_index
    _rnaseq_input = _taiji_input
    _rnaseq_output_dir = (<> "/RNASeq") . _taiji_output_dir
    _rnaseq_tmp_dir = _taiji_tmp_dir

instance SCRNASeqConfig TaijiConfig where
    _scrnaseq_batch_info = _taiji_batch_info
    _scrnaseq_genome_fasta = _taiji_genome
    _scrnaseq_star_index = _taiji_star_index
    _scrnaseq_annotation = fromMaybe
        (error "Please specificy the path of annotation file") . _taiji_annotation
    _scrnaseq_input = _taiji_input
    _scrnaseq_output_dir = (<> "/SCRNASeq") . _taiji_output_dir
    _scrnaseq_tmp_dir = _taiji_tmp_dir
    _scrnaseq_cell_barcode_length = _scrna_cell_barcode_length . _taiji_scrna_options
    _scrnaseq_molecular_barcode_length = _scrna_umi_length . _taiji_scrna_options
    _scrnaseq_doublet_score_cutoff = _scrna_doublet_score_cutoff . _taiji_scrna_options
    _scrnaseq_cluster_resolution_list = _scrna_cluster_resolution_list . _taiji_scrna_options
    _scrnaseq_cluster_resolution = _scrna_cluster_resolution . _taiji_scrna_options
    _scrnaseq_cluster_optimizer = _scrna_cluster_optimizer . _taiji_scrna_options

-- Construct workflow
build "wf" [t| SciFlow TaijiConfig |] $ do
    Core.builder
    SingleCell.builder

    namespace "RNA" RNASeq.builder
    namespace "ATAC" ATACSeq.builder
    [ "ATAC_Get_TFBS", "ATAC_Get_Peak", "HiC_Read_Input"
        , "RNA_Make_Expr_Table", "ATAC_Make_Expr_Table" ] ~> "Create_Linkage_Prep"
    ["Create_Linkage", "RNA_Make_Expr_Table"] ~> "Compute_Ranks_Prep"

    namespace "SCATAC" SCATACSeq.builder
    namespace "SCRNA" SCRNASeq.builder
    ["SCATAC_Find_TFBS", "SCATAC_Call_Peaks", "SCATAC_Gene_Acc"] ~>
        "Compute_Ranks_SC_Prep"

    Fusion.builder

getCoordConfig :: String -> Int -> FilePath -> IO RemoteConfig
getCoordConfig ip port fl = do
    config <- getDefaultRemoteConfig ["remote", "--ip", ip, "--port", show port]
    settings <- decodeFileThrow fl :: IO (M.HashMap String Value)
    return config
        { _remote_parameters = str <$> M.lookup "submit_params" settings
        , _submission_cmd = str $ M.lookupDefault "sbatch" "submit_command" settings
        , _cpu_format = str $ M.lookupDefault "--ntasks-per-node=%d" "submit_cpu_format" settings
        , _memory_format = str $ M.lookupDefault "--mem=%d000" "submit_memory_format" settings
        }
  where
    str (String x) = T.unpack x
    str _ = error "Expecting string"


main :: IO ()
main = defaultMain header descr commands wf
  where
    header = printf "Taiji-v%s" $ showVersion version
    descr = "Multi-omics bioinformatics analysis pipline. " <>
        "For more details, see: https://taiji-pipeline.github.io/"
    commands =
        [ runParser getCoordConfig
        , deleteParser
        , showParser
        , viewParser
        , remoteParser (Proxy :: Proxy Remote) ]
