Taiji -- multi-omics bioinformatics pipeline
============================================

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/Taiji-pipeline/Taiji?include_prereleases)
![Continuous integration](https://github.com/Taiji-pipeline/Taiji/workflows/Continuous%20integration/badge.svg)
![GitHub All Releases](https://img.shields.io/github/downloads/Taiji-pipeline/Taiji/total)

Taiji is an integrative analysis pipeline for analyzing bulk/single-cell ATAC-seq and RNA-seq data.
Please go to this [website](https://taiji-pipeline.github.io/) for documentation and tutorials. 

- Joint analysis of ATAC-seq, RNA-seq and Hi-C datasets.
- Integrate multiple single cell datasets, scale to more than 1 million cells!

The design philosophy of the Taiji pipeline is focused on:

- Correctness: We only include reliable algorithms and
    make every effort to ensure the implementations are bug-free.
- Performance: We code algorithms from scratch when necessary to
    ensure the pipeline can scale to large datasets (thousands of samples at least).
- Convinence: Most analyses have multipe entry points, e.g., Fastq, Bam or Bed.
    The execution of the pipeline requires only a single command.

We achieve these at the expense of customization. This will be improved in the future.

Installation
------------

Pre-built binaries are available for macOS and Linux system:

- `taiji-CentOS-x86_64`: for Red Hat Enterprise Linux derivatives.

- `taiji-Ubuntu-x86_64`: for Debian linux derivatives.

- `taiji-macOS-XX-XX`: for macOS.

Example:

```
curl -L https://github.com/Taiji-pipeline/Taiji/releases/latest/download/taiji-CentOS-x86_64 -o taiji
chmod +x taiji
./taiji --help
```

If you have used **Taiji** in your research, please consider citing the following paper:

K. Zhang, M. Wang, Y. Zhao, W. Wang,
Taiji: System-level identification of key transcription factors reveals
transcriptional waves in mouse embryonic development.
*Sci. Adv.* **5**, eaav3262 (2019).
