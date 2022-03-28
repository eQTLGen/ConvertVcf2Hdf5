#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="ConvertVcf2Hdf5"

#These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# Define paths
nextflow_path=../../tools/  # Path to Nextflow executable, no need to adjust if folder structure is same as recommended in cookbook.

genopath=[Folder with input genotype files in .vcf.gz format]
outputpath=../output/ # Path to output folder, no need to adjust if the folder structure is same as recommended in cookbook.
cohortname=[your study name CohortName_ExpressionPlatformName]

# Optional argument
# --samplelist '[Path to file with sample IDs]'

NXF_VER=20.10.0 ${nextflow_path}/nextflow run ConvertVcf2Hdf5.nf \
--vcf ${genopath} \
--outdir ${outputpath} \
--cohort_name ${cohortname} \
-profile slurm,singularity \
-resume
