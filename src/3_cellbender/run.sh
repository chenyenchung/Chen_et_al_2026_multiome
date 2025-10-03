#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --gres=gpu:1
#SBATCH --out=log/cellbender_%j.out

module purge
module load r/gcc/4.5.0
base_path=${1}

name=outs/raw_feature_bc_matrix.h5

cd ${base_path}

UTIL_PREFIX="/scratch/ycc520/thesis/src/cellbender/bin"
Rscript ${UTIL_PREFIX}/extractGEX.R --h5 "${name}"

IMG_PREFIX="/scratch/cgsb/desplan/File_exchange/Yen_ref"
singularity exec --nv \
  --overlay ${IMG_PREFIX}/overlay/cellbender.ext3:ro \
  ${IMG_PREFIX}/image/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif \
  bash -c "source /ext3/bin/enable_mamba.sh; micromamba run -n cellbender cellbender remove-background --cuda --input outs/gex --output outs/gex_cellbender.h5"
