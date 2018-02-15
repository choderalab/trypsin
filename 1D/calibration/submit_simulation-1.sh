#!/bin/bash 
#BSUB -J Imatinib-solvent-1-GBAOAB
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -q gpuqueue
#BSUB -R select[gpu_model0=='GeForceGTX1080Ti']
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
source activate protons-production

# Launch my program.
module load cuda
python run_simulation-1.py

# Submit the next job into the queue
bsub < submit_simulation-2.sh
