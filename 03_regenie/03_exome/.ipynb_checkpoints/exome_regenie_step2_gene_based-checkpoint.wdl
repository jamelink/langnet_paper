version 1.1

task exome_regenie_step2_gene_based {

    input {

        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File annotation_file
        File set_list
        File? exclude_sets
        Boolean firth = true
        File mask_definitions
        File phenotype_file
        File covariate_file
        File step1_predictions 
        File step1_predictions_list
        File extract_sets
        String gene_based_tests
        String build_mask
        String categorical_covariates
        String outname
        Float max_alternate_allele_frequency
        Int mac_threshold
        Int min_mac
        Boolean write_mask = false
        String? test_type
        String? aaf_bins
        String? skat_params

    }

    String plink_bed_basename = basename("${plink_bed_file}", ".bed")

    command <<<

    set -x -e -o pipefail

    # Move input files to a new directory
    mkdir /input_files

    mv ~{plink_bed_file} /input_files/
    mv ~{plink_bim_file} /input_files/
    mv ~{plink_fam_file} /input_files/
    mv ~{step1_predictions} /input_files/
    tar xfzC /input_files/step_1_sent_exome.tar.gz /input_files/

    regenie \
    --step 2 \
    --bed /input_files/~{plink_bed_basename} \
    --pred ~{step1_predictions_list} \
    --covarFile ~{covariate_file} \
    --catCovarList ~{categorical_covariates} \
    --phenoFile ~{phenotype_file} \
    --build-mask ~{build_mask} \
    --write-mask-snplist \
    --check-burden-files \
    --vc-tests ~{gene_based_tests} \
    --vc-maxAAF ~{max_alternate_allele_frequency} \
    --vc-MACthr ~{mac_threshold} \
    --minMAC ~{min_mac} \
    --anno-file ~{annotation_file} \
    --set-list ~{set_list} \
    --extract-sets ~{extract_sets} \
    --mask-def ~{mask_definitions} \
    --apply-rint \
    --gz \
    ~{if defined(exclude_sets) then "--exclude-sets ${exclude_sets} \\" else "\\"}
    ~{if firth then "--firth --approx \\" else "--spa \\"}
    ~{if defined(test_type) then "--test ${test_type} \\" else "\\"}
    ~{if defined(aaf_bins) then "--aaf-bins ${aaf_bins} \\" else "\\"}
    ~{if defined(skat_params) then "--skat-params ${skat_params} \\" else "\\"}
    ~{if write_mask then "--write-mask \\" else "\\"}
    --verbose \
    --out ~{outname}

    >>>
    
    output {

        Array[File] regenie_output = glob("~{outname}*.regenie.gz")
        Array[File] regenie_log = glob("~{outname}*.log")
        Array[File]? masks_output = glob("~{outname}_masks*")

    }

    runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/regenie_3.2.1.tar.gz"
    
    }


}

