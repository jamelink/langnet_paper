chr=$1
#ancestry=$2
set_no=$2

if [ -z ${chr} ]; then
echo "Please provide a chromosome number (1-22, X or Y) when running this script."
exit 1
fi

if [ -z ${set_no} ]; then
echo "Please provide a set number when running this script."
exit 1
fi

dx_indir="AICHA_rs_conn:/exome_filtering/filtered_plink_bgen"
jobfile_prefix="exome_regenie_gene_based_additive_sum_AICHA"
dx_outdir="AICHA_rs_conn:/association_analyses/output_gene_based/"
firth="false"
gene_based_tests="skat,skato"
build_mask="sum"
outname_prefix="AICHA_rs_conn_gene_based_additive_sum"
max_alternate_allele_frequency=0.01
mac_threshold=0
min_mac=1
write_mask="false"
test_type=""
aaf_bins=0.01
skat_params="1,25"

# Make a folder where jobs are stored
mkdir -p ./jobfiles

# Set different instances for chromosomes
#if [[ ${chr} == 2 ]] || [[ ${chr} == 19 ]]; then
#instance_type="mem3_ssd1_v2_x8"
#else
instance_type="mem1_ssd1_v2_x4"
#fi

echo -n "dx run \"AICHA_rs_conn:/workflows/exome_regenie_step2_gene_based/exome_regenie_step2_gene_based\" \\
--destination=\"${dx_outdir}/c${chr}/\" \\
--cost-limit=4 \\
--tag=\"regenie_step2_gene_based\" \\
--tag=\"c${chr}\" \\
--tag=\"additive_sum\" \\
--name=\"regenie_step2_gene_based_c${chr}_set_${set_no}\" \\
--priority=\"normal\" \\
--instance-type=\"${instance_type}\" \\
--yes \\
--brief \\
-iplink_bed_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.bed\" \\
-iplink_bim_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.bim\" \\
-iplink_fam_file=\"${dx_indir}/c${chr}/AICHA_rs_conn_ukb23157_c${chr}_v1_variant_filter.fam\" \\
-iannotation_file=\"AICHA_rs_conn:/exome_filtering/annotations/AICHA_rs_conn_c${chr}_annotation_file.txt\" \\
-iset_list=\"AICHA_rs_conn:/exome_filtering/annotations/AICHA_rs_conn_c${chr}_set_list.txt\" \\
-iextract_sets=\"AICHA_rs_conn:/exome_filtering/annotations/AICHA_rs_conn_c${chr}_extract_set_${set_no}.txt\" \\
-ifirth=${firth} \\
-imask_definitions=\"AICHA_rs_conn:/exome_filtering/annotations/handedness_mask_definitions.txt\" \\
-ibuild_mask=\"${build_mask}\" \\
-iphenotype_file=\"AICHA_rs_conn:/phenotype/sent_edges_exome.tsv\" \\
-icovariate_file=\"AICHA_rs_conn:/phenotype/covars_pc10_exome.tsv\" \\
-icategorical_covariates=\"sex,geno_array_dummy,site_dummy_11025,site_dummy_11026,site_dummy_11027,exome_batch\"  \\
-istep1_predictions=\"AICHA_rs_conn:/association_analyses/step1_predictions/step_1_sent_exome.tar.gz\" \\
-istep1_predictions_list=\"AICHA_rs_conn:/association_analyses/step1_predictions/reg_st1_exome_pred.list\" \\
-imax_alternate_allele_frequency=${max_alternate_allele_frequency} \\
-imin_mac=${min_mac} \\
-imac_threshold=${mac_threshold} \\
-igene_based_tests=\"${gene_based_tests}\" \\
-iwrite_mask=${write_mask} \\
-ioutname=\"${outname_prefix}_c${chr}_set_${set_no}\" \\
-iaaf_bins=\"${aaf_bins}\" \\
-iskat_params=\"${skat_params}\"

" > ./jobfiles/${jobfile_prefix}_c${chr}_set_${set_no}.sh

awk '{$1=$1;print}' ./jobfiles/${jobfile_prefix}_c${chr}_set_${set_no}.sh > tmp && mv tmp ./jobfiles/${jobfile_prefix}_c${chr}_set_${set_no}.sh

bash ./jobfiles/${jobfile_prefix}_c${chr}_set_${set_no}.sh