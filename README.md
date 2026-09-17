# Domain-resolved multiome profiling of the L3 OPC neurons

Code for processing paired RNA and chromatin-accessibility data, annotating cell populations, testing differential expression and accessibility, and computing corrected local Simpson indices (cLSI). The study and demonstration dataset is [**GSE344177**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344177). Analyses are implemented as separate workflows under `src/`.

## System requirements

**Target platform:** Linux HPC with Slurm

| Component | Version | Use |
| --- | --- | --- |
| Nextflow | 26.04.0 | Workflow execution; later workflows use typed processes and workflow outputs |
| OpenJDK | 21.0.3 | Nextflow runtime recorded in execution logs |
| R / renv | 4.5.1 / 1.1.5 | R analyses and package restoration |
| Bioconductor | 3.22 in `renv.lock` | Bioconductor package release |
| Seurat / SeuratObject / Signac | 5.4.0 / 5.3.0 / 1.16.0 | Multiome objects and analysis |
| edgeR / Matrix | 4.8.2 / 1.7-4 | Pseudobulk tests and sparse matrices |
| R.utils / data.table / ggplot2 | 2.13.0 / 1.18.0 / 4.0.1 | Script arguments, tables, and figures |
| RANN / ggrastr / dunn.test / tidyr | 2.6.2 / 1.0.2 / 1.3.7 / 1.3.2 | cLSI and statistical plots |
| karyoploteR / GenomicRanges / rtracklayer | 1.36.0 / 1.62.1 / 1.70.1 | Genomic interval plots |
| GenomeInfoDb / TxDb.Dmelanogaster.UCSC.dm6.ensGene | 1.46.2 / 3.12.0 | Genome annotations |
| patchwork / RColorBrewer | 1.3.2 / 1.1-3 | Figure assembly and palettes |
| Cell Ranger ARC | 2.0.2 | Multiome alignment |
| CellBender | 0.3.2 | Background removal; packaged in the GHCR image below |
| BWA-MEM2 / FastQC / Cutadapt | 2.2.1 / 0.11.9 / 4.9 | Whole-genome sequencing preparation |
| GATK / Picard | 4.3.0.0 / 2.27.5 | Variant calling |
| SAMtools / BCFtools | 1.14 and 1.20 / 1.14 and 1.19 | Versions vary by process |
| BEDTools / BEDOPS / VCFtools | 2.29.2 / 2.4.39 / 0.1.16 | Module versions in preprocessing |
| deepTools | 3.5.6 | Coverage tracks |
| Python | 2.7, 3.12, or 3.14 in separate process environments | MACS2 and metadata processing |
| Go | 1.20 | Building scbamop |

The full R dependency list and package sources are in [`renv.lock`](renv.lock). R 4.5.1 environment is expected.

Slurm, Git, environment modules, Micromamba, and Singularity are required. `conda-forge` version as of 2026-09-16 of MEME Suite, BEDTools, SAMtools, MACS2, pandas, SciPy, Scanpy, and Muon are used in workflows. CellBender and popscle use external environments or container images without complete build specifications in this repository. Exact versions must be recorded from the environments used for the submission.

The published CellBender image combines the CUDA 12.8.1/cuDNN 9.8.0 Ubuntu 24.04.2 base with an embedded ext3 overlay containing the `cellbender` Micromamba environment.

**Hardware:** CPU stages request up to 24 cores and 48 GB RAM per task. Figure 2 tasks request 1 core and 4–32 GB RAM each. Concurrent tasks require additional aggregate resources. CellBender requests an NVIDIA GPU. GPU model, GPU memory requirement, and peak disk usage have not been benchmarked. The full pipeline is configured for HPC; execution on a standard laptop has not been validated.

## Installation

1. Obtain the repository and run subsequent commands from its root. Initialize the pinned submodules. For HTTPS access without GitHub SSH keys:

   ```bash
   git -c url.https://github.com/.insteadOf=git@github.com: submodule update --init --recursive
   export PROJECT_ROOT="$PWD"
   mkdir -p log result int work
   ```

2. Load the site installations of Java 21, Nextflow 26.04.0, R 4.5.1, and Micromamba. Load Singularity for container-based stages. Module names are site-specific; on the development cluster the R and Nextflow modules are:

   ```bash
   module load r/4.5.1
   module load nextflow/26.04.0
   java -version
   nextflow -version
   Rscript --version
   micromamba --version
   ```

3. Restore the R environment; do not regenerate the lockfile:

   ```bash
   Rscript -e 'renv::restore(prompt = FALSE)'
   Rscript -e 'renv::status()'
   ```

   The project `.Rprofile` activates renv. Scripts launched from Nextflow task directories load `R.utils` before loading the project environment. Make the restored library visible to those processes:

   ```bash
   export R_LIBS_USER="$(Rscript -e 'cat(renv::paths$library())')"
   ```

4. For annotation and coverage export, build scbamop and expose it on `PATH`:

   ```bash
   (cd extern/scbamop && go build -o build/scbamop .)
   export PATH="$PROJECT_ROOT/extern/scbamop/build:$PATH"
   ```

   Pull the all-in-one CellBender SIF from GHCR (approximately 38 GB):

   ```bash
   mkdir -p "$PROJECT_ROOT/static"
   singularity pull "$PROJECT_ROOT/static/cellbender.sif" \
     oras://ghcr.io/chenyenchung/cellbender:v20260915
   ```

   The ext3 overlay is embedded in the SIF; no separate `--overlay` argument is needed. The SIF file's SHA-256 is `633e7f76814a8d6b3096994f03cad1c261112e016d1bba39bbab5a62d047bbb7`. Verify the installed environment with:

   ```bash
   singularity exec --cleanenv "$PROJECT_ROOT/static/cellbender.sif" \
     bash -c 'source /ext3/bin/enable_mamba.sh; micromamba run -n cellbender cellbender --version'
   ```

   Expected output: `0.3.2`.

5. Configure the scheduler in a local `hpc.config`, for example:

   ```groovy
   process.executor = 'slurm'
   process.queue = 'YOUR_PARTITION'
   process.clusterOptions = '--account=YOUR_ACCOUNT'
   ```

   Replace the example partition and account with site values. Preserve the GPU requests in GPU processes when adapting scheduler settings. The repository configuration enables Conda with Micromamba; process environments are created on first use. Compute nodes need access to the repository, restored R library, reference files, containers, and work directory.

**Expected installation time:** allow approximately 1–3 hours on a normal Linux machine with compilers installed and a working internet connection for R restoration and CPU dependencies. This is an unmeasured planning estimate and excludes data transfer, container provisioning, and scheduler delays. Installation and first-run timings must be measured before submission.

## Demo

### Data availability

Use the study data deposited under [GEO accession GSE344177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344177). Accession contents and download availability were not verified when preparing this README. No GSE344177 test subset, processed object, or expected numerical output is tracked in Git. The development checkout's `data` directory points to external project storage. A fresh clone therefore requires data staging before the dataset demonstration can run.

The smallest dataset demonstration below starts from the study's processed `annotated.rds` Seurat object, produced by the annotation stage. Obtain that object with its associated study metadata, or generate it using the upstream stages described below. Raw GEO downloads alone cannot be passed directly to this command.

### Run differential analysis

On a compute node with **1 CPU and 16 GB RAM**, after installation:

```bash
export PROJECT_ROOT="$PWD"
export DEMO_OBJECT=/absolute/path/to/GSE344177/annotated.rds
mkdir -p result/demo
(
  cd result/demo
  /usr/bin/time -p Rscript "$PROJECT_ROOT/src/7_analyses/bin/run_de.R" \
    --proot "$PROJECT_ROOT" \
    --obj "$DEMO_OBJECT" \
    --edger_helper "$PROJECT_ROOT/src/7_analyses/bin/edger_ql_helpers.R"
)
```

Expected files under `result/demo/`:

```text
de_identity/RNA.csv
de_identity/ATAC.csv
de_spatial/RNA.csv
de_spatial/ATAC.csv
```

Each CSV contains feature-level edgeR results, including log fold changes, P values, FDR, feature identifiers, contrasts, and cell-class labels. Identity tables also contain `contrast_type`. Identity analyses compare each level with the remaining levels; spatial tables contain pairwise comparisons. Reference row counts and numerical checksums are not supplied.

**Expected runtime:** unmeasured for GSE344177. The corresponding Nextflow process requests 30 minutes, 1 CPU, and 16 GB RAM; this is a scheduler limit, not a measured runtime. Record elapsed time from `/usr/bin/time` on the staged dataset.

An independent check using included synthetic labels is available:

```bash
Rscript src/common/test_clsi_helpers.R
```

Expected message: `All cLSI helper tests passed.` This check should take seconds (unmeasured estimate); it does not exercise the GSE344177 workflow.

## Instructions for use

### Entry points and input preparation

Run workflows from the repository root. Supply absolute input paths. Overriding `--root` does not rewrite other default paths or paths embedded in scripts. Inspect the parameters and process commands of each selected workflow before launch.

| Stage | Entry point | Inputs and purpose |
| --- | --- | --- |
| WGS preparation | `src/1_preparation/align_wgs.nf`, `call_var.nf` | Paired FASTQs, genome FASTA and indexes, variant resources; alignment and genotype calling |
| Multiome alignment | `src/1_preparation/align_arc.nf` | `data/stf_*.csv` Cell Ranger ARC library sheets and `--ref` reference directory |
| Demultiplexing | `src/2_preprocessing/demuxlet.nf` | Aligned libraries, DGRP genotypes, reporter-line variants, line list, popscle image |
| Background removal | `src/3_cellbender/run.sh` | Library directory containing `outs/raw_feature_bc_matrix.h5` |
| Object preparation | `src/4_obj_prep/main.nf` | Aligned libraries, demultiplexing results, barcode filters, genome annotations, marker and reference assets |
| Annotation | `src/5_cluster_annotation/main.nf` | Clustered object, insertion files, BAMs, barcodes, and genome annotations |
| Figure 2 analysis | `src/7_analyses/fig2.nf` | Annotated object, reference genome, motifs, consensus peaks, and external reference objects |
| ATAC coverage export | `src/7_analyses/export_bw.nf` | Cluster CSV, per-library BAM directories, and scbamop executable |

For a new dataset, create Cell Ranger ARC library sheets with the required `fastqs,sample,library_type` fields. Adapt library lists, barcode conventions, genotype resources, reference paths, and annotation mappings in the relevant stages. Use a consistent genome assembly and chromosome naming throughout; the study defaults use Drosophila BDGP6/dm6 and Ensembl release 88.

Before submitting `src/3_cellbender/run.sh` with `sbatch`, change its `src/cellbender/bin` path to the actual `src/3_cellbender/bin` path. Replace its Singularity invocation with the following, using the absolute path to the downloaded SIF and running from the library directory after GEX extraction:

```bash
singularity exec --nv /absolute/path/to/thesis/static/cellbender.sif \
  bash -c 'source /ext3/bin/enable_mamba.sh; micromamba run -n cellbender cellbender remove-background --cuda --input outs/gex --output outs/gex_cellbender.h5'
```

### Reproduce quantitative analyses and figures

The Figure 2 workflow computes differential expression/accessibility, motif enrichment, cLSI comparisons, QC, and figure panels. It also requires external Özel reference objects and an origin mapping; these assets are not tracked in Git. With all inputs staged, run:

```bash
export PROJECT_ROOT="$PWD"
export STUDY_INPUTS=/absolute/path/to/staged/study_inputs
nextflow run src/7_analyses/fig2.nf -c hpc.config \
  --root "$PROJECT_ROOT" \
  --obj "$STUDY_INPUTS/annotated.rds" \
  --fasta "$STUDY_INPUTS/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa" \
  --gtf "$STUDY_INPUTS/Drosophila_melanogaster.BDGP6.88.gtf" \
  --all_peaks "$STUDY_INPUTS/consensus_peaks.bed" \
  --motif "$STUDY_INPUTS/CisBPDrosophilaALL.meme" \
  --origin_map "$STUDY_INPUTS/origin_map.csv" \
  --ozel_obj_dir "$STUDY_INPUTS/ozel_2020_obj" \
  --edger_helper "$PROJECT_ROOT/src/7_analyses/bin/edger_ql_helpers.R" \
  --clsi_helper "$PROJECT_ROOT/src/common/clsi_helpers.R" \
  -with-report log/fig2-report.html \
  -with-trace log/fig2-trace.tsv \
  -with-timeline log/fig2-timeline.html
```

Outputs are copied beneath `result/`, grouped by the workflow's published output names. They include the four differential-analysis tables; XSTREME output; `ISS_raw.csv`, `clsi_spatial_vs_temporal.csv`, and `clsi_spatial_over_class.csv`; and PDF panels for DEG/DAR counts, genomic locations, cLSI comparisons, QC, and UMAP. Exact directory mappings are defined in the `publish` and `output` blocks of `fig2.nf`.

The differential analysis aggregates counts by cell class, identity, and library, applies TMM normalization and expression filtering, and fits edgeR models with biological replicate adjustment. Spatial and temporal tests use a log2 fold-change threshold of 0.5; Notch tests use 1. Figure count plots select FDR < 0.05 and exclude Notch comparisons. Preserve these settings and the original biological replicate assignments when reproducing the study.

Coverage export accepts `--cluster_csv`, `--bam`, and `--scbamop`; BAMs must follow `<bam-root>/<library>/outs/atac_possorted_bam.bam`.

