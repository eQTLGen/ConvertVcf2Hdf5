# Genotype conversion pipeline

This is Nextflow pipeline for converting imputed .vcf file folder into hdf5 format. Genotypes in hdf5 format can be used as an input for HASE tool.

## Usage instructions

### Requirements for the system

- Have access to HPC with multiple cores.

- Have Bash >=3.2 installed.

- Have Java 8 installed.

- Have SLURM scheduler managing the jobs in the HPC.

- Have Singularity (preferred) and/or conda installed.

### Requirements for the input data

- Input .vcf folder should contain 22 bgzipped .vcf.gz files.

- It should not contain other file types (e.g. .vcf.gz.tbi or README).

- Currently the pipeline assumes that there is genotype probabilities in the second field for each individual (GT:**GP**).

### Setup the pipeline

1. Make a folder named `genotype_conversion` and get in there.

2. Make subfolder `tools`, go into this folder and have Nextflow executable installed: https://www.nextflow.io/docs/latest/getstarted.html#installation. After that, step back up into `genotype_conversion`.

3. Get the genotype conversion pipeline from here: https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5

You can either clone it by using git (if available in HPC):

`git clone https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5`

Or just download this from the gitlab download link and unzip into `genotype_conversion` folder.

4. Make folder `help_files`.

5. Get or make the file with SNP IDs to filter in and put it into `help_files` folder. SNP format has to be: chr[chrnum]\_[hg19pos]\_[refallele]\_[altallele]. For eQTLGen phase 2 analyses, get file **SNPList_MAF_0005.txt**: this includes SNPs with MAF>0.005 in 1000G p3v5 ALL reference panel. 

6. Make the tab-separated file which lists all the samples to include into hdf5 file and put it into `help_files`. For eQTLGen phase 2 analyses, this list should contain only genotype IDs for those samples which are included into final eQTL mapping.

7. If your system does not support Singularity and you use conda instead, download Genotype-IO from here and put it into folder `/genotype_conversion/ConvertVcf2Hdf5/bin/`. If you use Singularity, no action is needed as this tool is already present in the container.

8. Make folder where you want to save your genotype data in hdf5 format. This can be in the `genotype_conversion` folder or elsewhere.

### Running the conversion command

Go to folder `ConvertVcf2Hdf5` and modify the SLURM script template `submit_genotype_conversion_template.sh` with your full input paths.

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
outputpath=[Folder where to write genotype data in h5 format]
snplist=[path to the SNP ID list]
studyname=[your study nameCohortName_ExpressionPlatformName]

# Optional argument for the command
# --samplelist '[Path to file with sample IDs]'

NXF_VER=20.10.0 ${nextflow_path}/nextflow run ConvertVcf2Hdf5.nf \
--inputpath ${genopath} \
--outputpath ${outputpath} \
--snplist ${snplist} \
--studyname ${studyname} \
-profile slurm,singularity \
-with-report GenotypeConversionReport.html \
-resume
```

You can save the modified script version to informative name, e.g. `submit_genotype_conversion_[**CohortName_PlatformName**].sh`.

Then submit the job `sbatch submit_genotype_conversion_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.

### Monitoring and cleaning up

- Monitor the `slurm-***.out` log file and check if all steps finish without error.

	Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.

- If pipeline crashes and you figure out how to fix this, you can just resubmit the same script after fixes. Nextflow does not rerun unaffected steps and continues only from the steps affected by the issue. 

- When the work has finished, check `GenotypeConversionReport.html` for potential errors or warnings.

- For some analyses, `work` directory can become quite large. So, when the pipeline is successfully finished, you can delete it in order to save disk space of your HPC. But this means that the pipeline will restart from scratch, if you ever need to rerun it.

### Results of the pipeline

After successful completion of the pipeline, there should be hdf5 genotype file structure in your output folder, with folders named `individuals`, `probes`, `genotype` and `SNPQC`. For eQTLGen phase 2 analyses, this folder is one of the inputs for PerCohortPreparations pipeline.

### Benchmark

Here is the estimate, how much time the conversion is expected to take.

- Initial number of samples in .vcf genotype data: \~10,000

- Number of samples to filter into .h5 genotype format: \~1,000

- Infrastructure: University HPC with \~150 compute nodes

- Dependency management: Singularity 

- Time to run the pipeline (without wall times): \~4.5h

- CPU hours: \~80h

- Final size of `work` subdirectory: \~500GB

## Acknowledgements

This pipeline utilizes HASE (https://github.com/roshchupkin/hase) and some of its helper scripts originally developed by:

Gennady V. Roscupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands) 

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands). 

### Changes from the original HASE repo

Robert Warmerdam (Department of Genetics, University Medical Center Groningen, University of Groningen, Groningen, Netherlands) modified the original HASE and fixed some bugs.

Urmo Võsa (Institute of Genomics, University of Tartu, Tartu, Estonia) incorporated it into Nextflow pipeline and customized some parts of the code.

**Changes:**

- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the `--intercept` option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Lines 355-357 of hase.py were commented out because this caused pipeline to crash when >1 datasets were added.
- Line 355 of /hdgwas/data.py were changed `self.chunk_size=10000 --> self.chunk_size=20000`.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep the tool lightweight.

### Citation

Original method paper for HASE framework:

[Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)](https://www.nature.com/articles/srep36076)

### Contacts

For this Nextflow pipeline: urmo.vosa@gmail.com

For the method of HASE, find contacts from here: https://github.com/roshchupkin/hase
