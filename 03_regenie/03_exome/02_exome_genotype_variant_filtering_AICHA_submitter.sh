chr=$1

if [ -z ${chr} ]; then
echo "Please provide a chromosome number (1-22, X or Y) when running this script."
exit 1
fi

vcf_blocks=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/helper_files/pvcf_blocks.txt
dx_bulk_dir="AICHA_rs_conn:/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"
dx_outdir="AICHA_rs_conn:/exome_filtering/output"

#Extract the chromosome from the vcf blocks file
blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} )

# Make a folder where jobs are stored
mkdir -p ./jobfiles/genotype_variant_filtering

for block in ${blocks}; do

    echo -n "dx run \"AICHA_rs_conn:/workflows/exome_genotype_variant_filtering/exome_genotype_variant_filtering\" \\
    --destination=\"${dx_outdir}/c${chr}/\" \\
    --cost-limit=7 \\
    --tag=\"genotype_variant_filtering\" \\
    --tag=\"c${chr}\" \\
    --tag=\"b${block}\" \\
    --name=\"genotype_variant_filtering_c${chr}_b${block}\" \\
    --priority=\"normal\" \\
    --yes \\
    --brief \\
    -istage-common.subject_list=\"AICHA_rs_conn:file-GQ48v8jJZzvPx9QFypJGY9pZ\" \\
    -istage-common.target_regions=\"AICHA_rs_conn:file-GFV6B5jJxJ1JJZjZPKvKFf9P\" \\
    -istage-common.project_name=\"AICHA_rs_conn\" \\
    -istage-common.DP_snv=7 \\
    -istage-common.DP_indel=10 \\
    -istage-common.GQ=20 \\
    -istage-common.average_GQ=35 \\
    -istage-common.min_MAC=1 \\
    -istage-common.genotype_missingness=0.10 \\
    -istage-common.min_AB_snv=0.15 \\
    -istage-common.min_AB_indel=0.2 \\
    -istage-common.vcf_file=\"${dx_bulk_dir}/ukb23157_c${chr}_b${block}_v1.vcf.gz\" \\
    -istage-common.vcf_index=\"${dx_bulk_dir}/ukb23157_c${chr}_b${block}_v1.vcf.gz.tbi\"
    " > ./jobfiles/genotype_variant_filtering/exome_genotype_variant_filtering_submit_c${chr}_b${block}.sh

    awk '{$1=$1;print}' ./jobfiles/genotype_variant_filtering/exome_genotype_variant_filtering_submit_c${chr}_b${block}.sh > tmp && mv tmp ./jobfiles/genotype_variant_filtering/exome_genotype_variant_filtering_submit_c${chr}_b${block}.sh

    bash ./jobfiles/genotype_variant_filtering/exome_genotype_variant_filtering_submit_c${chr}_b${block}.sh

done
