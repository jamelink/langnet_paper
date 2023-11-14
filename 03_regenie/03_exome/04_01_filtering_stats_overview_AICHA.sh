vcf_blocks=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/helper_files/pvcf_blocks.txt
chr=$1
workdir=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/wes/
pvcfsdir=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/

mkdir -p ${workdir}/genotype_variant_filtering/c${chr}
cd ${workdir}/genotype_variant_filtering/c${chr}

dx cd /exome_filtering/output/c${chr}/

dx download -a "AICHA_rs_conn_ukb23157_c${chr}_b*_v1_pre_filtering_summary.txt"
dx download -a "AICHA_rs_conn_ukb23157_c${chr}_b*_v1_variant_filter_summary.txt"

dx cd /exome_filtering/filtered_plink_bgen/c${chr}

dx download -a "AICHA_rs_conn_ukb23157_c${chr}_b*_v1_variant_filter_conversion_summary.txt"

blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} )

for block in ${blocks}; do

printf "Collecting statistics for chromosome ${chr}, block ${block}..."

echo ${chr} >> ${workdir}/genotype_variant_filtering/c${chr}/chromosome.txt
echo ${block} >> ${workdir}/genotype_variant_filtering/c${chr}/block.txt

awk '$1 == 3 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_pre_filtering_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/variants_pre.txt

awk '$1 == 3 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/avg_gq.txt
awk '$1 == 4 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/missing.txt
awk '$1 == 5 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/mac.txt
awk '$1 == 6 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/ab.txt
awk '$1 == 7 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/total_removed.txt
awk '$1 == 8 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/ts_pre.txt
awk '$1 == 9 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/tv_pre.txt
awk '$1 == 10 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/ts_post.txt
awk '$1 == 11 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/tv_post.txt
awk '$1 == 12 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/total_remaining.txt

awk '$1 == 3 {print $3}' ${workdir}/genotype_variant_filtering/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_b${block}_v1_variant_filter_conversion_summary.txt >> ${workdir}/genotype_variant_filtering/c${chr}/multiallelic_removed.txt

module load vcftools/0.1.17
module load samtools/1.9



gzip -cd ${pvcfsdir}/sites_only_pvcfs/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.vcf.gz | grep -v "#" | wc -l >> ${workdir}/genotype_variant_filtering/c${chr}/total_all.txt

gzip -cd ${pvcfsdir}/sites_only_pvcfs/c${chr}/ukb23157_c${chr}_b${block}_v1_target_site_only.vcf.gz | grep -v "#" | wc -l >> ${workdir}/genotype_variant_filtering/c${chr}/total_target_regions.txt

gzip -cd ${pvcfsdir}/sites_only_pvcfs/c${chr}/ukb23157_c${chr}_b${block}_v1_target_site_only.vcf.gz | grep -v "#" | awk '$7 == "MONOALLELIC"' | wc -l >> ${workdir}/genotype_variant_filtering/c${chr}/monoallelic.txt

printf "Done\n"

done

paste ${workdir}/genotype_variant_filtering/c${chr}/chromosome.txt ${workdir}/genotype_variant_filtering/c${chr}/block.txt ${workdir}/genotype_variant_filtering/c${chr}/total_all.txt ${workdir}/genotype_variant_filtering/c${chr}/total_target_regions.txt ${workdir}/genotype_variant_filtering/c${chr}/monoallelic.txt ${workdir}/genotype_variant_filtering/c${chr}/variants_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/avg_gq.txt ${workdir}/genotype_variant_filtering/c${chr}/missing.txt ${workdir}/genotype_variant_filtering/c${chr}/mac.txt ${workdir}/genotype_variant_filtering/c${chr}/ab.txt ${workdir}/genotype_variant_filtering/c${chr}/total_removed.txt ${workdir}/genotype_variant_filtering/c${chr}/total_remaining.txt ${workdir}/genotype_variant_filtering/c${chr}/ts_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/tv_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/ts_post.txt ${workdir}/genotype_variant_filtering/c${chr}/tv_post.txt ${workdir}/genotype_variant_filtering/c${chr}/multiallelic_removed.txt > ${workdir}/genotype_variant_filtering/c${chr}/c${chr}_overview_filtering_statistics_per_block.txt

rm ${workdir}/genotype_variant_filtering/c${chr}/chromosome.txt ${workdir}/genotype_variant_filtering/c${chr}/block.txt ${workdir}/genotype_variant_filtering/c${chr}/total_all.txt ${workdir}/genotype_variant_filtering/c${chr}/total_target_regions.txt ${workdir}/genotype_variant_filtering/c${chr}/monoallelic.txt ${workdir}/genotype_variant_filtering/c${chr}/variants_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/avg_gq.txt ${workdir}/genotype_variant_filtering/c${chr}/missing.txt ${workdir}/genotype_variant_filtering/c${chr}/mac.txt ${workdir}/genotype_variant_filtering/c${chr}/ab.txt ${workdir}/genotype_variant_filtering/c${chr}/total_removed.txt ${workdir}/genotype_variant_filtering/c${chr}/total_remaining.txt ${workdir}/genotype_variant_filtering/c${chr}/ts_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/tv_pre.txt ${workdir}/genotype_variant_filtering/c${chr}/ts_post.txt ${workdir}/genotype_variant_filtering/c${chr}/tv_post.txt ${workdir}/genotype_variant_filtering/c${chr}/multiallelic_removed.txt

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$13/$14,$15,$16,$15/$16,$17}' ${workdir}/genotype_variant_filtering/c${chr}/c${chr}_overview_filtering_statistics_per_block.txt > ${workdir}/genotype_variant_filtering/c${chr}/tmp && mv ${workdir}/genotype_variant_filtering/c${chr}/tmp ${workdir}/genotype_variant_filtering/c${chr}/c${chr}_overview_filtering_statistics_per_block.txt