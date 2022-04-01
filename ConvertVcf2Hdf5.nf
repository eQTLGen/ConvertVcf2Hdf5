#!/usr/bin/env nextflow

def helpmessage() {

log.info"""
ConvertVcfToHdf5 v${workflow.manifest.version}"
===========================================================
Pipeline for converting bgzipped .vcf.gz files into hdf5 format.

This pipeline makes extensive use of the help scripts from R.Roschupkin (https://github.com/roshchupkin/hase/) but, in order to speed those up, parallelizes conversion task in the HPC cluster.

Usage:

nextflow run ConvertVcfToHdf5.nf \
--vcf '/inputfolder/' \
--cohort_name 'StudyName' \
--outdir '/outputfolder/'

Mandatory arguments:
--vcf             Path to input vcf folder. It must contain bgzipped .vcf.gz files.
--outdir          Path to output folder where genotype data is written in .hdf5 format.
--cohort_name     Name of the study.

Optional arguments:
--VcfOutput     If directory is specified, outputs filtered vcf files into this folder.
""".stripIndent()

}

// Default parameters
params.vcf = ''
params.outdir = ''
params.cohort_name = ''
params.VcfOutput = ''

//Show parameter values
log.info """================================================================
ConvertVcfToHdf5 v${workflow.manifest.version}"
================================================================"""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Input genotype directory']                 = params.vcf
summary['Output directory']                         = params.outdir
summary['Cohort name']                              = params.cohort_name
summary['Output directory for vcf']                 = params.VcfOutput

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.vcf}/*chr${chr}.filtered.vcf.gz")) }
  .ifEmpty { exit 1, "Input .vcf.gz files not found!" } 
  .into { chr_vcf_pairs_fixnames; chr_vcf_pairs_count_samples }

process CountSamples {

    input:
      tuple chr, file(vcf) from chr_vcf_pairs_count_samples

    output:
       path "sample_list.txt" into SampleList

    when:
      chr==22

    script:
      """
      # Count the number of samples in the vcf file
      bcftools query -l ${vcf} > sample_list.txt
      """
}

NumberOfSamples = SampleList.splitCsv().count().get()
log.info "Number of detected samples: " + NumberOfSamples

process FixVariantNames {

    tag {"FixVariantNames_$chr"}

    input:
      tuple chr, file(vcf) from chr_vcf_pairs_fixnames

    output:
      tuple chr, file("${chr}_FixedSnpNames.vcf.gz") into vcf_fixed_snp_names_ch
      file("${chr}_FixedSnpNames.vcf.gz") into VcfToMakeIndAndProbe 

    script:
      """
      bcftools view --max-alleles 2 ${vcf} \
      | awk 'BEGIN {FS="\\t"; OFS="\\t"}; {if (substr(\$0,0,1) !~ "#"){\$3="chr"  \$1  ":"  \$2  "_"  \$4  "_" \$5;}print \$0;}' \
      | bgzip -c > ${chr}_FixedSnpNames.vcf.gz
      """
}

vcf_fixed_snp_names_ch.into{VcfToRemoveInfo; VcfToTabix}

process RemoveInfoField {

    tag {"RemoveInfo_$chr"}

    input:
      tuple chr, file(vcf) from VcfToRemoveInfo

    output:
      tuple chr, file() into VcfToChunkVcf

    script:
      """
      bcftools annotate -x INFO ${vcf} | bgzip -c > ${chr}_FixedSnpNamesInfoRemoved.vcf.gz
      """
}


process MakeIndAndProbe {

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
      path InputFiles from VcfToMakeIndAndProbe.collect()

    output:
      path "probes" 
      path "individuals"

    script:
      """
      mkdir -p InputPath
      mv ${InputFiles} InputPath/.

      bash $baseDir/bin/helperscripts/VCF2hdf5_SortedChr.sh InputPath . $baseDir/bin/hase/ ${params.cohort_name}

      # Remove unnecessary files
      rm -rf ./tmp_files
      rm -rf ./*.txt
      """
}

process ChunkVcf {

    tag {"ChunkVcf_$chr"}

    input:
      tuple chr, file(vcf) from VcfToChunkVcf

    output:
      path ("chr${chr}*") into dosagechunks

    script:
      """
      bcftools query -f "[\\t%DS]\\n" ${vcf} |\
      cut -f2- |\
      split -a 3 -d -l 25000 - chr${chr}_dosage
      """
}

process FixChunkSize {

    input:
      file input_dosage_chunks from dosagechunks.flatten().toSortedList()

    output:
      file ('*.txt') into DosageChunksToHdf5

    script:
      """
      Rscript --vanilla $baseDir/bin/helperscripts/ChunkFixer.R \$PWD 25000

      # Rename the chunks
      ind=0
      for i in `ls \$PWD | sort -V`
      do

      echo \${i}

      mv \${i} \${ind}_${params.cohort_name}.txt
      ind=\$((\$ind + 1))

      done
      """
}

process ConvertGenotypeToHdf5 {

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
      file (InpHdf5Chunks) from DosageChunksToHdf5.flatMap()

    output:
      file ('genotype/*.h5') 

    script:
      """
      echo ${InpHdf5Chunks}
      chunk_id=\$(echo ${InpHdf5Chunks} | sed -e "s/_.*//")

      mkdir output

      python $baseDir/bin/hase/tools/VCF2hdf5.py -flag chunk -id \${chunk_id} -data ${InpHdf5Chunks} -out output -study_name ${params.cohort_name}

      mkdir genotype
      mv output/genotype/*.h5 genotype/
      """
}

process TabixFilteredVcfInput {

    tag {"TabixFilteredVcfInput_$chr"}

    input:
      tuple chr, file(InputVcf) from VcfToTabix

    output:
      tuple chr, file(InputVcf), file("${InputVcf}.tbi") into (InputToSnpQc,VcfToOutput)

    script:
      """
      tabix -p vcf ${InputVcf}
      """
}

process CalculateSnpQcMetrics {

    tag {"CalculateSnpQcMetrics_$chr"}

    input:
      tuple chr, file(InputToSnpQc), file(InputToSnpQc_index) from InputToSnpQc

    output:
      file ("${chr}_statistics.vars") into SNPQC_files
    
    script:
      """
      java -Xmx8g -Xms8g -jar $baseDir/bin/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar \
        -i ${InputToSnpQc} \
        -I VCF \
        -o ${chr}_statistics

      bcftools query -f "%IMPUTED\\t%TYPED\\t%TYPED_ONLY\\n" ${InputToSnpQc} > ${chr}_imputation_info
      paste -d '\t' ${chr}_statistics.vars ${chr}_imputation_info > ${chr}_statistics.vars
      """
}

process CompressSnpQcFile {

    publishDir "${params.outdir}/SNPQC", mode: 'copy', overwrite: true

    input:
      file SnpQcReport from SNPQC_files.collectFile(name: "${params.cohort_name}_SNPQC.txt", keepHeader: true, sort: true)

    output:
      path("${params.cohort_name}_SNPQC.txt.gz")

    script:
      """
      gzip -f ${SnpQcReport}
      """
}

process OutputVcf {

    publishDir "${params.VcfOutput}", mode: 'copy', overwrite: true

    input:
      file(filtred_vcf) from VcfToOutput

    output:
      file(filtred_vcf) into OuputVcf

    when:
      params.VcfOutput

    script:
      """
      echo ${filtred_vcf}
      """

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
