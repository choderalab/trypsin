#!/bin/bash
#BSUB -J {0}
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -R select[gpu_model0=='GeForceGTX1080Ti']
#BSUB -gpu num=1:mode=exclusive_process:mps=no:
#BSUB -W  120:00
#BSUB -We 119:30
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

# quit on first error
set -e

# Change to working directory used for job submission
cd $LS_SUBCWD
export PATH="/home/rustenburg/miniconda3/bin:$PATH"

# Use the right conda environment
source activate trypsin

# Launch my program.
module load cuda/9.0
python {1} {2}
