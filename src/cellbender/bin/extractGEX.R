#!/usr/bin/env Rscript
library(Seurat)  
library(Matrix)
library(R.utils)

# Source: https://github.com/cellgeni/cellbender/blob/fdc46038eff406e4f346207ce014a817cbfd54ca/scripts/multiome.R

# this scripts goes through subfolders in data, extracts gene expression part of multiome 
# and saves it into gex subfolder as matrix.mtx.gz accompained by features.tsv.gz and barcodes.tsv.gz
# it should works with both h5 and mtx.gz input (but I didn't test it with the second)
# Pasha M

subsetGEX = function(d,out,type='Gene Expression'){
  print(d)
  m = Seurat::Read10X_h5(d,unique.features = FALSE,use.names = FALSE)[[type]]
  n = Seurat::Read10X_h5(d,unique.features = FALSE,use.names = TRUE)[[type]]
  features = data.frame(rownames(m),rownames(n),'Gene Expression')
  barcodes = colnames(m)
  
  if(!dir.exists(out)) dir.create(out)
    
  gz1 = gzfile(file.path(out,'/features.tsv.gz'), "w")
  write.table(features, gz1, quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  close(gz1)
      
  gz1 = gzfile(paste0(out,'/barcodes.tsv.gz'), "w")
  writeLines(barcodes,gz1)
  close(gz1)
            
  if(file.exists(paste0(out,'/matrix.mtx.gz'))) file.remove(paste0(out,'/matrix.mtx.gz'))
  writeMM(m,paste0(out,'/matrix.mtx')) # it doesn't work (or was extremely slow) with gzfile probably because of files sizes 
  system({paste0("gzip ",out,'/matrix.mtx')})
}

path = R.utils::commandArgs(asValues = TRUE)[["h5"]]

subsetGEX(path, paste0(dirname(path),'/gex'))
