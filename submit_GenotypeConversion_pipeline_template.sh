#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="ConvertVcf2Hdf5"

# These are needed modules in UT HPC to get Singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook
nextflow_path=../../tools/  # Path to Nextflow executable, no need to adjust if folder structure is same as recommended in cookbook.

cohortname=[name of your cohort]
genopath=../../2_Imputation/output   # Folder with input genotype files in .vcf.gz format
outputpath=../output/ # Path to output folder, no need to adjust if the folder structure is same as recommended in cookbook.

NXF_VER=21.10.6 ${nextflow_path}/nextflow run ConvertVcf2Hdf5.nf \
--vcf ${genopath} \
--outdir ${outputpath} \
--cohort_name ${cohortname} \
-profile slurm,singularity \
-resume
