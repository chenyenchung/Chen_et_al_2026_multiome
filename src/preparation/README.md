## General

Please refer to the `module` directives in the `nextflow` scripts and set
your environment up with the required dependencies.

## ./renv_setup.R

Setup R environment for reproducibility.

## ./align_wgs.nf

TODO: Give an accession number for WGS fastqs

Align whole genome sequencing results of heterozygous reporter lines
(crossed with reference genome Dm6 lines).

[bwa2-mem](https://github.com/bwa-mem2/bwa-mem2)v2.2.1 was installed
locally and required to reproduce this pipeline.

## ./call_var.nf

Call variants against Ensembl88 reference genome with bam files generated
by `align_wgs.nf`.

## ./align_arc.nf

Align multiome results.

[cellranger-arc](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest)
v2.0.2 was installed locally.

