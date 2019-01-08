#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=myMatlabTest
#SBATCH --mail-type=END
#SBATCH --mail-user=jmn425@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2018b

cd /scratch/$USER/LDA_NMF
cat hpccode.m | srun matlab -nodisplay

