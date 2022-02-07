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
nextflow_path=[full path to your Nextflow executable]

genopath=[Folder with input genotype files in .vcf.gz format]
outputpath=[Folder where to write genotype data in h5 format]
snplist=[path to the SNP ID list]
studyname=[your study nameCohortName_ExpressionPlatformName]

# Optional argument for the command
# --samplelist '[Path to file with sample IDs]'

NXF_VER=20.10.0 ${nextflow_path}/nextflow run ConvertVcf2Hdf5.nf \
--vcf ${genopath} \
--outdir ${outputpath} \
--snplist ${snplist} \
--studyname ${studyname} \
-profile slurm,singularity \
-with-report GenotypeConversionReport.html \
-resume
