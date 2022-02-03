
refdir=[folder with reference]
outdir=[output folder]

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' ${refdir}CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz > ${outdir}1Kg_30x_SNP_info.txt

for i in {2..22}
do
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' ${refdir}CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz >> ${outdir}1Kg_30x_SNP_info.txt
echo ${i}

done

