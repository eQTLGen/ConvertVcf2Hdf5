# Genotype conversion pipeline

This Nextflow pipeline can be used to convert imputed .vcf file folder into hdf5 format, usable by HASE tool.

## Requirements for the system

- Have access to HPC with multiple cores.

- Have bash installed.

- Have Java 8 installed.

- Have SLURM scheduler managing the jobs in the HPC.

- Have Singularity (preferred) or conda installed.

## Requirements for the input data

- Input .vcf folder should contain 22 .vcf or bgzipped .vcf.gz files.

- It should not contain other file types (e.g. .vcf.gz.tbi or README).

- Currently the pipeline assumes that there is genotype probabilities in the second field for each individual (GT:**GP**).

## Setup the pipeline

- Make a folder `genotype_conversion` and access there.

- Make folder `tools`, go into this folder and have Nextflow executable installed: https://www.nextflow.io/docs/latest/getstarted.html#installation. After that, go back to `genotype_conversion`.

- Get the genotype conversion pipeline from here: https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5

You can either clone it by using git (if available):

`git clone https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5`

Or just download this from the link and unzip to `genotype_conversion` folder.

- Make folder `help_files`.

- Get file **SNPList_MAF_0005.txt** from here and put it into `help_files` folder.

- Make the tab-separated file listing all the samples to include into output hdf5 file and put it into `help_files`.

## Running the conversion command

Go to folder `ConvertVcf2Hdf5` and modify the script template `submit_genotype_conversion_template.sh` with your inputs.

```
#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunHASEwNextflow"

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tools depend on HPC.
module load java-1.8.0_40
module load python/2.7.15/native
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=../tools

${nextflow_path}/nextflow run test_ConvertVcf2Hdf5.nf \
--inputpath '[**/path/to/your/input/genotype/vcf/folder/**]' \
--outputpath '[**/path/to/your/output/genotype/hdf5/folder/**]' \
--snplist '../help_files/SNPList_MAF_0005.txt' \
--samplelist '../help_files/[**FileWithYourSampleIds.txt**]' \
--studyname '[**CohortName_PlatrformName**]' \
-with-report GenotypeConversionReport.html \
-resume \
-profile [**singularity_profile**/**conda_profile**],cluster_slurm

```
You can save the modified script version to e.g. `submit_genotype_conversion_[**CohortName_PlatrformName**].sh`.

Then submit the job `sbatch submit_genotype_conversion_[**CohortName_PlatrformName**].sh`.

Monitor the `slurm-***.out` log file and check if all steps finish.

When work has finished, check `GenotypeConversionReport.html` to potential errors or warnings.



