#!/bin/bash
#SBATCH -J HaplaMiscarrRound3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16100M
#SBATCH --output=/ifs/data/research/CGM/Projects/Miscarr/1_Paper/Families/ErrorOutput/MiscarrFileConv.out
#SBATCH --error=/ifs/data/research/CGM/Projects/Miscarr/1_Paper/Families/ErrorOutput/MiscarrFileConv.err
#SBATCH --partition=defq

module load singularity

singularity exec /ifs/data/research/CGM/Containers/PGT/pgtR-1.0.simg Rscript /ifs/data/research/CGM/Projects/Miscarr/1_Paper/Fambuilder_2020.R 
