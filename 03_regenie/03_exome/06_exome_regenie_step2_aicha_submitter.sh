chr=2
range="2:24030937-24039490"
gene_name="WDCP"


if [ -z ${chr} ]; then
echo "Please provide a chromosome number (1-22, X or Y) when running this script."
exit 1
fi

#vcf_blocks=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/handedness/wes_data_450k/helper_files/pvcf_blocks.txt
dx_indir="AICHA_rs_conn:/exome_filtering/filtered_plink_bgen"
dx_outdir="AICHA_rs_conn:/association_analyses/output_sig_variants/"
#dx_helpers="AICHA_rs_conn:/exome_filtering/helper_files/sig_variant_ids"

#Extract the chromosome from the vcf blocks file

#blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} )
#last_block=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${vcf_blocks} | tail -1 )

# Make a folder where jobs are stored
mkdir -p ./jobfiles/regenie/

echo -n "dx run \"AICHA_rs_conn:/workflows/exome_regenie_step2/exome_regenie_step2_AICHA\" \\
-iplink_bed_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.bed\" \\
-iplink_bim_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.bim\" \\
-iplink_fam_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.fam\" \\
-iphenotype_file=\"AICHA_rs_conn:/phenotype/sent_edges_exome.tsv\" \\
-icovariate_file=\"AICHA_rs_conn:/phenotype/covars_pc10_exome.tsv\" \\
-istep1_predictions=\"AICHA_rs_conn:/association_analyses/step1_predictions/step_1_sent_exome.tar.gz\" \\
-ipredictions_file=\"AICHA_rs_conn:/association_analyses/step1_predictions/reg_st1_exome_pred.list\"  \\
-icategorical_covariates=\"sex,geno_array_dummy,site_dummy_11025,site_dummy_11026,site_dummy_11027,exome_batch\"  \\
-ioutname=\"sent_all_st2_gene_vars_${gene_name}\"  \\
-ivariant_range=\"${range}\" \\
-imin_mac=1 \\
--destination=\"${dx_outdir}/${gene_name}\" \\
--cost-limit=3 \\
--tag=\"regenie_st2\" \\
--tag=\"c${chr}\" \\
--name=\"regenie_st2_${gene_name}\" \\
--priority=\"normal\" \\
--instance-type=\"mem1_ssd1_v2_x4\" \\
--yes \\
--brief \\

" > ./jobfiles/regenie/exome_regenie_step2_aicha_${gene_name}.sh

awk '{$1=$1;print}' ./jobfiles/regenie/exome_regenie_step2_aicha_${gene_name}.sh > tmp && mv tmp ./jobfiles/regenie/exome_regenie_step2_aicha_${gene_name}.sh

bash ./jobfiles/regenie/exome_regenie_step2_aicha_${gene_name}.sh
