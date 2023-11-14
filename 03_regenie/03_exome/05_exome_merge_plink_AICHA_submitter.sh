chr=$1

if [ -z ${chr} ]; then
echo "Please provide a chromosome number (1-22, X or Y) when running this script."
exit 1
fi

vcf_blocks=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/helper_files/pvcf_blocks.txt
dx_indir="AICHA_rs_conn:/exome_filtering/filtered_plink_bgen"
dx_outdir="AICHA_rs_conn:/exome_filtering/filtered_plink_bgen"
project_name="AICHA_rs_conn"

#Extract the chromosome from the vcf blocks file
blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} )
last_block=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} | tail -1 )

# Make a folder where jobs are stored
mkdir -p ./jobfiles/merge_plink

echo -n "dx run \"${project_name}:/workflows/exome_merge_plink/exome_merge_plink_files\" \\
--destination=\"${dx_outdir}/c${chr}/\" \\
--cost-limit=5 \\
--tag=\"merge_plink\" \\
--tag=\"c${chr}\" \\
--name=\"merge_plink_c${chr}\" \\
--priority=\"normal\" \\
--yes \\
--brief \\
-istage-common.out_prefix=\"${project_name}_ukb23157_c${chr}_v1_variant_filter\" \\
" > ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

if [[ ${chr} == "X" ]]; then

    echo -n "-istage-common.split_x=true \\
    -istage-common.individuals_sex=\"AICHA_rs_conn:/phenotype/sex_file.txt\" \\
    " >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

else

    echo -n "-istage-common.split_x=false \\
    "  >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

fi


for block in ${blocks}; do

    if [ ${block} -eq ${last_block} ]; then

        echo "-istage-common.plink_bed_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.bed\" \\" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh
        echo "-istage-common.plink_bim_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.bim\" \\" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh
        echo "-istage-common.plink_fam_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.fam\"" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

    else

        echo "-istage-common.plink_bed_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.bed\" \\" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh
        echo "-istage-common.plink_bim_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.bim\" \\" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh
        echo "-istage-common.plink_fam_files=\"${dx_indir}/c${chr}/${project_name}_ukb23157_c${chr}_b${block}_v1_variant_filter.fam\" \\" >> ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

    fi

done

awk '{$1=$1;print}' ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh > tmp && mv tmp ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh

bash ./jobfiles/merge_plink/exome_merge_plink_c${chr}.sh
