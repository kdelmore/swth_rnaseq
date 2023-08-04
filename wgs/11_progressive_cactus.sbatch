#!/bin/bash
#SBATCH --job-name=align
#SBATCH --ntasks=48
#SBATCH --time=48:00:00
#SBATCH --mem=5000GB
#SBATCH --output=align.%j
#SBATCH --partition=big-mem

module load GCC/7.3.0-2.30  OpenMPI/3.1.1 Python/3.7.0

source /scratch/user/delmore/tools/cactus-bin-v2.6.4/cactus_env/bin/activate

export TMPDIR=/scratch/user/delmore/cnes/cactus_tmp

cactus \
        /scratch/user/delmore/cnes/jobstore \
        /scratch/user/delmore/cnes/seqfile \
        /scratch/user/delmore/cnes/test.hal \
        --batchSystem single_machine \
        --maxCores 48 \
        --maxMemory 500G \
        --restart
