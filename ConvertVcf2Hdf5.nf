#!/usr/bin/env nextflow

def helpmessage() {

log.info"""

ConvertVcfToHdf5Parallelized v${workflow.manifest.version}"
===========================================================
Pipeline for converting .vcf or bgzipped .vcf.gz files into hdf5 format.

This pipeline makes extensive use of the help scripts from R.Roschupkin (https://github.com/roshchupkin/hase/) but, in order to speed those up, parallelizes conversion task in the HPC cluster.

Usage:

nextflow run ConvertVcfToHdf5Parallelized.nf --inputpath '/inputfolder/' --outputpath '/outputfolder/' --studyname 'StudyName' --snplist [file with SNP list]

Mandatory arguments:
--inputpath     Path to input vcf folder. It must contain only .vcf files or bgzipped .vcf.gz files. Must not contain any other file types (i.e. tabix .tbi).
--outputpath    Path to output folder where genotype data is written in .hdf5 format.
--studyname     Name of the study.
--snplist       SNPs to include into analysis.
Optional arguments:
--samplelist    Samples to include into analysis.
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

//Show parameter values
log.info """================================================================
Genotype converter v${workflow.manifest.version}"
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

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

InpDir = Channel.fromPath(params.inputpath)
OutpDir = Channel.fromPath(params.outputpath)
Studyname = Channel.value(params.studyname)
SNPlist = file(params.snplist)
if (params.samplelist) {SampleList = file(params.samplelist)}

OrigVcf = Channel.fromPath(params.inputpath + "/*chr*")
OrigVcf.into{OrigVcfToFilterSamples; OrigVcfToRenameFile}

process FilterSamples {

    tag {FilterSamples}

    cpus 1
    memory '15 GB'
    time '12h'
    executor 'slurm'
    clusterOptions '--job-name=FilterSamples'

    input:
      path vcf from OrigVcfToFilterSamples
      file samplefile from SampleList

    output:
      env (chr) into (ChrFilterSamplesToFixVariantNames, ChrFilterSamplesToFilterVariants)
      file ('*_SamplesFiltered.vcf.gz') into VcfFilterSamplesSamplesFiltered

    when:
        params.samplelist

    shell:
        '''
        chr="$(echo !{vcf} |\
        sed -e "s/.*chr/chr/g" |\
        grep -oP "chr[0-9]{1,2}")"

        echo ${chr}

        bcftools view \
        -S !{samplefile} \
        --force-samples \
        !{vcf} | bgzip -c > ${chr}_SamplesFiltered.vcf.gz
        '''
}

process RenameFiles {

    tag {RenameFile}

    cpus 1
    memory '15 GB'
    time '12h'
    executor 'slurm'
    clusterOptions '--job-name=RenameFiles'

    input:
      path vcf from OrigVcfToRenameFile

    output:
      env (chr) into (ChrRenameFileToFixVariantNames, ChrRenameFileToFilterVariants)
      file ('*_SamplesFiltered.vcf.gz') into VcfRenameFilesSamplesFiltered

    when:
        !params.samplelist

    shell:
        '''
        chr="$(echo !{vcf} |\
        sed -e "s/.*chr/chr/g" |\
        grep -oP "chr[0-9]{1,2}")"

        echo ${chr}

        # Filter in subset of samples
        mv !{vcf} ${chr}_SamplesFiltered.vcf.gz
        '''
}

process FixVariantNames {

    tag {FixVariantNames}

    cpus 1
    memory '15 GB'
    time '12h'
    executor 'slurm'
    clusterOptions '--job-name=FixVariantNames'

    input:
      path vcf from VcfFilterSamplesSamplesFiltered.mix(VcfRenameFilesSamplesFiltered)

    output:
      file ('*_FixedSnpNames.vcf.gz') into VcfFixedSnpNames
      env (chr) into ChrToFilterVariants

 
    """
    chr="\$(echo ${vcf} |\
    sed -e "s/.*chr/chr/g" |\
    grep -oP "chr[0-9]{1,2}")"

    echo \${chr}

    bcftools view --max-alleles 2 ${vcf} \
    | awk 'BEGIN {FS="\\t"; OFS="\\t"}; {if (substr(\$0,0,1) !~ "#"){\$3="chr"  \$1  ":"  \$2  "_"  \$4  "_" \$5;}print \$0;}' \
    | bgzip -c > \${chr}_FixedSnpNames.vcf.gz
    """
}

process FilterVariants {
    tag {FilterVariants}

    cpus 1
    memory '10 GB'
    time '10h'
    executor 'slurm'
    clusterOptions '--job-name=FilterVariants'

    input:
      path vcf from VcfFixedSnpNames
      file InclusionList from SNPlist
      env chr from ChrToFilterVariants

    output:
      file('*_FixedSnpNamesFiltered0005.vcf.gz') into (VcfToChunkVcf, VcfToMakeIndAndProbe, VcfToCountFiles, VcfToTabix, FixedNamesFilteredToHdfInput)

    """    
    vcftools --gzvcf ${vcf} \
    --snps ${InclusionList} \
    --recode \
    --stdout | bgzip -c \
    > \${chr}_FixedSnpNamesFiltered0005.vcf.gz
    """
}

CountFiles = VcfToCountFiles.collect().size()
InputFilesToIndProbe = VcfToMakeIndAndProbe.collect()

process MakeIndAndProbe {

    tag {MakeIndAndProbe}

    cpus 1
    memory '15 GB'
    time '12h'
    executor 'slurm'
    clusterOptions '--job-name=MakeIndAndProbe'

    input:
      path InputFiles from InputFilesToIndProbe
      val NameOfStudy from Studyname
      val CountFiles from CountFiles

    output:
      path OutputPath into IndProbeOutput

    when:
      CountFiles == 22

    """
    InpPath="InputPath"
    OutPath="OutputPath"

    mkdir -p \${InpPath}
    mkdir -p \${OutPath}

    mv ${InputFiles} \${InpPath}/.

    bash $baseDir/bin/helperscripts/VCF2hdf5_SortedChr.sh \${InpPath} \${OutPath} $baseDir/bin/hase/ ${NameOfStudy}

    # Remove unnecessary files
    rm -rf \${OutPath}/tmp_files
    rm -rf \${OutPath}/*.txt

    # Make separate folder for SNP QC file
    mkdir \${OutPath}/SNPQC
    """
}

process ChunkVcf {

    tag {ChunkVcf}

    cpus 1
    memory '5 GB'
    time '2h'
    executor 'slurm'
    clusterOptions '--job-name=ChunkVcf'

    input:
      path vcf from VcfToChunkVcf

    output:
      path ('chr*') into rawchunks

    """
    chr="\$(echo ${vcf} |\
    sed -e "s/.*chr/chr/g" |\
    grep -oP "chr[0-9]{1,2}")"

    echo \${chr}

    zcat -f ${vcf} |\
    awk 'BEGIN{FS="\t"}/^[^#]/{print }' |\
    cut -f10- |\
    split -a 3 -d -l 25000 - \${chr}_
    """
}

process CalculateDosage {

    tag {CalculateDosage}

    cpus 1
    memory '5 GB'
    time '45min'
    executor 'slurm'
    clusterOptions '--job-name=CalculateDosage'

    input:
      file input_dosage_chunks from rawchunks.flatten()

    output:
      file ('chr*_dosage') into dosagechunks

    """
    echo ${input_dosage_chunks}
    cat ${input_dosage_chunks} | awk 'BEGIN{FS="\t"}{R=""; for (i=1; i<=NF; i++){split(\$i, a, ":"); split(a[2], b, ","); if(i==1){R=b[2] + 2*b[3]}else{R=R"\t"b[2] + 2*b[3]}}; print R}' > ${input_dosage_chunks}_dosage
    """
}

process FixChunkSize {

    tag {FixChunkSize}

    cpus 4
    memory '25 GB'
    time '6h'
    executor 'slurm'
    clusterOptions '--job-name=FixChunkSize'

    input:
      file input_dosage_chunks from dosagechunks.toSortedList()
      val (NameOfStudy) from Studyname

    output:
      file ('*.txt') into DosageChunksToHdf5

    """
    Rscript --vanilla $baseDir/bin/helperscripts/ChunkFixer.R \$PWD 25000

    # Rename the chunks
    ind=0
    for i in `ls \$PWD | sort -V`
    do

    echo \${i}

    mv \${i} \${ind}_${NameOfStudy}.txt
    ind=\$((\$ind + 1))

    done
    """
}

process ConvertGenotypeToHdf5 {

  tag {ConvertGenotypeToHdf5}

    cpus 1
    memory '2 GB'
    time '10min'
    executor 'slurm'
    clusterOptions '--job-name=ConvertGenToHdf5'

    input:
      file (InpHdf5Chunks) from DosageChunksToHdf5.flatMap()
      val (NameOfStudy) from Studyname

    output:
      file ('output/genotype/*.h5') into Hdf5ToCollect

    """
    echo ${InpHdf5Chunks}
    chunk_id=\$(echo ${InpHdf5Chunks})

    mkdir output

    python $baseDir/bin/hase/tools/VCF2hdf5.py -flag chunk -id \${chunk_id} -data ${InpHdf5Chunks} -out output -study_name ${NameOfStudy}
    
    """
}

process CollectHdf5Chunks {

  tag {CollectHdf5Chunks}

    cpus 1
    memory '2 GB'
    time '10min'
    executor 'slurm'
    clusterOptions '--job-name=CollectHdf5Chunks'

    input:
      file (Hdf5Files) from Hdf5ToCollect.flatMap()
      path IndProbeFolder from IndProbeOutput

    output:
      path IndProbeFolder into ProbeAndIndToSnpQc

    """
    mv ${Hdf5Files} ${IndProbeFolder}/genotype/.
    """
}

process TabixFilteredVcfInput {

    tag {TabixFilteredVcfInput}

    cpus 1
    memory '10 GB'
    time '2h'
    executor 'slurm'
    clusterOptions '--job-name=TabixFilteredVcfInput'

    input:
      file InputVcf from VcfToTabix

    output:
      tuple file(InputVcf), file("*_FixedSnpNamesFiltered0005.vcf.gz.tbi") into InputToSnpQc
      file "*_FixedSnpNamesFiltered0005.vcf.gz.tbi" into TbiIndexFileToCount

    """
    tabix -p vcf ${InputVcf}
    """
}

process CalculateSnpQcMetrics {

    tag {CalculateSnpQcMetrics}

    cpus 1
    memory '10 GB'
    time '10h'
    executor 'slurm'
    clusterOptions '--job-name=CalculateSnpQcMetrics'

    input:
      path InputToSnpQc from InputToSnpQc

    output:
      file ('*.vars') into (SNPQC_files, SNPQC_files2)
    
    script:
    if(workflow.containerEngine == 'singularity'){
    """
    chr="\$(echo ${InputToSnpQc} |\
    sed -e "s/.*chr/chr/g" |\
    grep -oP "chr[0-9]{1,2}")"

    echo chr\${chr}

    java -Xmx10g -Xms10g -jar /usr/bin/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar \
    -i ${InputToSnpQc}/chr\${chr}_FixedSnpNamesFiltered0005.vcf.gz \
    -I VCF \
    -o chr\${chr}_statistics
    """
    } else {
    """
    chr="\$(echo ${InputToSnpQc} |\
    sed -e "s/.*chr/chr/g" |\
    grep -oP "chr[0-9]{1,2}")"

    echo chr\${chr}

    java -Xmx10g -Xms10g -jar $baseDir/bin/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar \
    -i ${InputToSnpQc}/chr\${chr}_FixedSnpNamesFiltered0005.vcf.gz \
    -I VCF \
    -o chr\${chr}_statistics
    """
    }
}


HelpChannel3 = SNPQC_files2.collect()
CountFiles3 = HelpChannel3.size()

process CompressSnpQcFile {

    tag {CompressSnpQcFile}

    cpus 1
    memory '5 GB'
    time '1h'
    executor 'slurm'
    clusterOptions '--job-name=CompressSnpQcFile'

    publishDir "${params.outputpath}", mode: 'copy', overwrite: true

    input:
      file SnpQcReport from SNPQC_files.flatten().collectFile(name: 'SNPQC.txt', keepHeader: true, sort: true)
      val CountFiles3 from CountFiles3
      path outputpath from OutpDir
      path ProbeAndIndToSnpQc from ProbeAndIndToSnpQc
      val (NameOfStudy) from Studyname

    output:
      path("*") into ProbeAndIndToWriteOut

    when:
      CountFiles3 == 22

    """
    gzip -f ${SnpQcReport}

    mv SNPQC.txt.gz ${ProbeAndIndToSnpQc}/SNPQC/${NameOfStudy}_SNPQC.txt.gz

    mv ${ProbeAndIndToSnpQc}/* .
    rm -r ${ProbeAndIndToSnpQc}

    """
}
