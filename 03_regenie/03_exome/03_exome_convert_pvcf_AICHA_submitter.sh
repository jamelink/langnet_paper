chr=$1

if [ -z ${chr} ]; then
echo "Please provide a chromosome number (1-22, X or Y) when running this script."
exit 1
fi

vcf_blocks=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/helper_files/pvcf_blocks.txt
dx_indir="AICHA_rs_conn:/exome_filtering/output"
dx_outdir="AICHA_rs_conn:/exome_filtering/filtered_plink_bgen"
project_name="AICHA_rs_conn"

#Extract the chromosome from the vcf blocks file
blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} )

# Make a folder where jobs are stored
mkdir -p ./jobfiles/convert_plink/

for block in ${blocks}; do

    echo -n "dx run \"${project_name}:/workflows/exome_convert_pvcf/exome_convert_pvcf_to_plink\" \\
    --destination=\"${dx_outdir}/c${chr}/\" \\
    --cost-limit=5 \\
    --tag=\"convert_pvcf_to_plink\" \\
    --tag=\"c${chr}\" \\
    --tag=\"b${block}\" \\
    --name=\"convert_pvcf_to_plink_c${chr}_b${block}\" \\
    --priority=\"normal\" \\
    --yes \\
    --brief \\
    -istage-common.vcf_file=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.vcf.gz\" \\
    -istage-common.vcf_index=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.vcf.gz.tbi\" \\
    " > ./jobfiles/convert_plink/exome_convert_pvcf_to_plink_submit_c${chr}_b${block}.sh

    awk '{$1=$1;print}' ./jobfiles/convert_plink/exome_convert_pvcf_to_plink_submit_c${chr}_b${block}.sh > tmp && mv tmp ./jobfiles/convert_plink/exome_convert_pvcf_to_plink_submit_c${chr}_b${block}.sh

     bash ./jobfiles/convert_plink/exome_convert_pvcf_to_plink_submit_c${chr}_b${block}.sh

done

