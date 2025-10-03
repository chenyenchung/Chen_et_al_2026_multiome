library(BSgenomeForge)
library(AnnotationForge)
library(rtracklayer)

# BSGenome
BSgenomeForge::fastaTo2bit(
  "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa",
  "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.2bit"
)

forgeBSgenomeDataPkgFromTwobitFile(
  "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.2bit",
  organism = "Drosophila melanogaster", provider = "BDGP",
  pkg_maintainer = "Yen-Chung Chen <yenchung.chen@nyu.edu>",
  circ_seqs = character(0),
  genome = "dm6",
  destdir = "/scratch/cgsb/desplan/File_exchange/Yen_ref/BSGenome_ensembl_88/"
)

# OrgDb
gtf <- import("/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf")
gtfdf <- as.data.frame(subset(gtf, type == "gene"))
fsym <- gtfdf[, c("gene_id", "gene_name")]
colnames(fsym) <- c("GID", "SYMBOL")
fchr <- gtfdf[, c("gene_id", "seqnames")]
colnames(fchr) <- c("GID", "CHROMOSOME")
fsym$ENSEMBL <- fsym$GID

makeOrgPackage(
  gene_info=fsym, chromosome=fchr,
  version = "1.0",
  maintainer = "Yen-Chung Chen <yenchung.chen@nyu.edu>",
  author = "Yen-Chung Chen <yenchung.chen@nyu.edu>",
  outputDir = "/scratch/cgsb/desplan/File_exchange/Yen_ref/BSGenome_ensembl_88/",
  tax_id = "7227",
  genus = "Drosophila",
  species = "melanogaster"
)

renv::install("/scratch/cgsb/desplan/File_exchange/Yen_ref/BSGenome_ensembl_88/org.Dmelanogaster.eg.db/")
renv::install("/scratch/cgsb/desplan/File_exchange/Yen_ref/BSGenome_ensembl_88/BSgenome.Dmelanogaster.BDGP.dm6/")
renv::snapshot()
