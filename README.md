# Genotype conversion pipeline

This is Nextflow pipeline for converting imputed .vcf file folder into hdf5 format. Genotypes in hdf5 format can be used as an input for HASE tool.

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java 8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone TBA`

Or just download this from the gitlab/github download link and unzip.

### Input files

#### Required

- Input .vcf folder should contain 22 bgzipped .vcf.gz files.
- Pipeline assumes that there is dosage (DS) field available in the imputed `.vcf` file. This is the case for the output of eQTLGen imputation pipeline.

#### Optional

- Flag `--samplelist` enables to specify file with the list of sample IDs to include into conversion. If you follow eQTLGen phase II cookbook, this is not crucial: only genotype samples which are in the genotype-to-expression file are included to imputation and `.hdf5` conversion steps. 

### Running the genotype conversion command

Go to folder `ConvertVcf2Hdf5` and modify the Slurm script template `submit_GenotypeConversion_pipeline_template.sh` with your input paths. This is an example template using Slurm scheduler.

```
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
outputpath=[Folder where to write genotype data in .hdf5 format]
cohort_name=[your study nam CohortName_ExpressionPlatformName]

# Optional argument for the command
# --samplelist '[Path to file with sample IDs]'

NXF_VER=21.10.6 ${nextflow_path}/nextflow run ConvertVcf2Hdf5.nf \
--vcf ${genopath} \
--outdir ${outputpath} \
--cohort_name ${cohort_name} \
-profile slurm,singularity \
-resume
```

You can save the modified script version to informative name, e.g. `submit_GenotypeConversion_pipeline_[**CohortName_PlatformName**].sh`.

Then submit the job `sbatch submit_GenotypeConversion_pipeline_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.

### Monitoring and debugging

- Monitoring:
  - Monitor the `slurm-***.out` log file and check if all the steps finish without error. Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.
  - Use `squeue -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file  is automatically written to your output folder `pipeline_info` subfolder, for potential errors or warnings. E.g. `output/pipeline_info/DataQcReport.html`.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from `work` folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.
  

### Output

After successful completion of the pipeline, there should be hdf5 genotype file structure in your output folder, with folders named `individuals`, `probes`, `genotype` and `SNPQC`. Additionally, there is folder `pipeline_info` with Nextflow runtime logs that can be used for debugging. For eQTLGen phase 2 analyses, this folder is one of the inputs for `PerCohortPreparations` pipeline.

## Acknowledgements

This pipeline utilizes HASE (https://github.com/roshchupkin/hase) and some of its helper scripts originally developed by:

Gennady V. Roscupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands) 

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands). 

### Changes from the original HASE repo

Robert Warmerdam (Department of Genetics, University Medical Center Groningen, University of Groningen, Groningen, Netherlands) modified the original HASE and fixed some bugs.

Urmo Võsa (Institute of Genomics, University of Tartu, Tartu, Estonia) incorporated it into Nextflow pipeline and customized/supplanted some parts of the code.

**Changes:**

- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the `--intercept` option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Setting /hdgwas/data.py were changed `self.chunk_size=10000 --> self.chunk_size=20000`.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep the tool lightweight.

### Citation

Original method paper for HASE framework:

[Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)](https://www.nature.com/articles/srep36076)

### Contacts

For this Nextflow pipeline: urmo.vosa at gmail.com

For the method of HASE, find contacts from [original HASE repo](https://github.com/roshchupkin/hase)
