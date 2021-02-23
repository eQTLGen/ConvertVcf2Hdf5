#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunHASEwNextflow"

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load java-1.8.0_40
module load python/2.7.15/native
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=../tools/

NXF_VER=20.10.0 ${nextflow_path}/nextflow run test_ConvertVcf2Hdf5.nf \
--inputpath '[Folder with genotype files in .vcf format]' \
--outputpath '[Folder where to write genotype data in h5 format]' \
--snplist '../help_files/SNPList_MAF_0005.txt' \
--samplelist '../help_files/[File with sample IDs]' \
--studyname '[CohortName_ExpressionPlatformName]' \
-profile [**singularity_profile**/**conda_profile**],cluster_slurm] \
-with-report GenotypeConversionReport.html \
-resume
