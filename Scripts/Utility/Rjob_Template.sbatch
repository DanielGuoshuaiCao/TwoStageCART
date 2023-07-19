#!/bin/bash
#SBATCH --job-name=IP_scenic
#SBATCH -o %x.%A.o
#SBATCH -e %x.%A.e
#SBATCH --time=20:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --mail-type=END
#SBATCH --mail-user=guoshuaicao@uchicago.edu

module load R/4.0.0
module load rstudio

cd /home/guoshuaicao/daniel/workspace/P3_CART_Traj/Scripts
Rscript Aphe_Regulons_cluster_version.R


