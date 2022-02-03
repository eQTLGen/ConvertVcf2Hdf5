#!/usr/bin/env nextflow

def helpmessage() {

log.info"""
ConvertVcfToHdf5 v${workflow.manifest.version}"
===========================================================
Pipeline for converting .vcf or bgzipped .vcf.gz files into hdf5 format.

This pipeline makes extensive use of the help scripts from R.Roschupkin (https://github.com/roshchupkin/hase/) but, in order to speed those up, parallelizes conversion task in the HPC cluster.

Usage:

nextflow run ConvertVcfToHdf5.nf --inputpath '/inputfolder/' --outputpath '/outputfolder/' --studyname 'StudyName' --snplist [file with SNP list]

Mandatory arguments:
--inputpath     Path to input vcf folder. It must contain bgzipped .vcf.gz files.
--outputpath    Path to output folder where genotype data is written in .hdf5 format.
--studyname     Name of the study.
--snplist       SNPs to include into analysis.
Optional arguments:
--samplelist    Samples to include into analysis.
--VcfOutput     If directory is specified, outputs filtered vcf files into this folder.
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.inputpath = ''
params.outputpath = ''
params.studyname = ''
params.samplelist = ''
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
summary['Input genotype directory']                 = params.inputpath
summary['Output directory']                         = params.outputpath
summary['SNP list']                                 = params.snplist
summary['Sample list']                              = params.samplelist
summary['Study name']                               = params.studyname
summary['Output directory for vcf']                 = params.VcfOutput

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.inputpath}/*chr${chr}_*.vcf.gz")) }
  .into { chr_vcf_pairs_rename_ch; chr_vcf_pairs_filter_ch }

snp_list_ch = Channel.fromPath(params.snplist).collect()

if (params.samplelist) { 
    sample_list_ch = Channel.fromPath(params.samplelist).collect()
    NumberOfSamples = Channel.fromPath(params.samplelist, checkIfExists: true).splitCsv().count().get()
    log.info "Number of detected samples: " + NumberOfSamples
  } else {
    NumberOfSamples = 5000
}

process FilterSamples {

    tag {"FilterSamples_$chr"}

    input:
      tuple chr, file(vcf) from chr_vcf_pairs_filter_ch
      file samplefile from sample_list_ch

    output:
      tuple chr, file("${chr}_SamplesFiltered.vcf.gz") into filtered_vcf_samples_ch

    when:
      params.samplelist

    script:
      """
      bcftools view \
      -S ${samplefile} \
      --force-samples \
      ${vcf} | bgzip -c > ${chr}_SamplesFiltered.vcf.gz
      """
}

process RenameFiles {

    tag {"RenameFile_$chr"}

    input:
      tuple chr, file(vcf) from chr_vcf_pairs_rename_ch

    output:
      tuple chr, file("${chr}_SamplesFiltered.vcf.gz") into renamed_vcf_samples_ch

    when:
      !params.samplelist

    script:
      """
      # Filter in subset of samples
      mv ${vcf} ${chr}_SamplesFiltered.vcf.gz
      """
}

process FixVariantNames {

    tag {"FixVariantNames_$chr"}

    input:
      tuple chr, file(vcf) from filtered_vcf_samples_ch.mix(renamed_vcf_samples_ch)

    output:
      tuple chr, file("${chr}_FixedSnpNames.vcf.gz") into vcf_fixed_snp_names_ch

    script:
      """
      bcftools view --max-alleles 2 ${vcf} \
      | awk 'BEGIN {FS="\\t"; OFS="\\t"}; {if (substr(\$0,0,1) !~ "#"){\$3="chr"  \$1  ":"  \$2  "_"  \$4  "_" \$5;}print \$0;}' \
      | bgzip -c > ${chr}_FixedSnpNames.vcf.gz
      """
}

process FilterVariants {

    tag {"FilterVariants_$chr"}

    input:
      tuple chr, file(vcf) from vcf_fixed_snp_names_ch
      file inclusion_snp_list from snp_list_ch

    output:
      tuple chr, file("${chr}_FixedSnpNamesFiltered0005.vcf.gz") into (VcfToChunkVcf, VcfToTabix)
      file("${chr}_FixedSnpNamesFiltered0005.vcf.gz") into VcfToMakeIndAndProbe
    
    script:
      """    
      vcftools --gzvcf ${vcf} \
      --snps ${inclusion_snp_list} \
      --recode \
      --stdout | bgzip -c \
      > ${chr}_FixedSnpNamesFiltered0005.vcf.gz
      """
}

process MakeIndAndProbe {

    publishDir "${params.outputpath}/", mode: 'copy', overwrite: true

    input:
      path InputFiles from VcfToMakeIndAndProbe.collect()

    output:
      path "probes" 
      path "individuals"

    script:
      """
      mkdir -p InputPath
      mv ${InputFiles} InputPath/.

      bash $baseDir/bin/helperscripts/VCF2hdf5_SortedChr.sh InputPath . $baseDir/bin/hase/ ${params.studyname}

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
      path ('chr*') into rawchunks

    script:
      """
      bcftools query -f "[\\t%DS]\\n" ${vcf} |\
      cut -f2- |\
      split -a 3 -d -l 25000 - chr${chr}_
      """
}

process CalculateDosage {

    input:
      file input_dosage_chunks from rawchunks.flatten()

    output:
      file ('chr*_dosage') into dosagechunks

    script:
      """
      cat ${input_dosage_chunks} > ${input_dosage_chunks}_dosage
      """
}

process FixChunkSize {

    input:
      file input_dosage_chunks from dosagechunks.toSortedList()

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

      mv \${i} \${ind}_${params.studyname}.txt
      ind=\$((\$ind + 1))

      done
      """
}

process ConvertGenotypeToHdf5 {

    publishDir "${params.outputpath}/", mode: 'copy', overwrite: true

    input:
      file (InpHdf5Chunks) from DosageChunksToHdf5.flatMap()

    output:
      file ('genotype/*.h5') 

    script:
      """
      echo ${InpHdf5Chunks}
      chunk_id=\$(echo ${InpHdf5Chunks} | sed -e "s/_.*//")

      mkdir output

      python $baseDir/bin/hase/tools/VCF2hdf5.py -flag chunk -id \${chunk_id} -data ${InpHdf5Chunks} -out output -study_name ${params.studyname}

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

    publishDir "${params.outputpath}/SNPQC", mode: 'copy', overwrite: true

    input:
      file SnpQcReport from SNPQC_files.collectFile(name: "${params.studyname}_SNPQC.txt", keepHeader: true, sort: true)

    output:
      path("${params.studyname}_SNPQC.txt.gz")

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
