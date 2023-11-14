version 1.1

task exome_regenie_step2_AICHA {

    input {

        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File phenotype_file
        File covariate_file
        File predictions_file
        File step1_predictions
        String variant_range
        String categorical_covariates
        String outname
        String min_mac

    }

    String plink_bed_basename = basename("${plink_bed_file}", ".bed")

    command <<<

    # Move input files to a new directory
    mkdir /input_files

    mv ~{plink_bed_file} /input_files/
    mv ~{plink_bim_file} /input_files/
    mv ~{plink_fam_file} /input_files/
    mv ~{step1_predictions} /input_files/
    tar xfzC /input_files/step_1_sent_exome.tar.gz /input_files/
    ls /input_files/*
    
    regenie \
    --step 2 \
    --bsize 200 \
    --bed /input_files/~{plink_bed_basename} \
    --pred ~{predictions_file} \
    --covarFile ~{covariate_file} \
    --catCovarList ~{categorical_covariates} \
    --phenoFile ~{phenotype_file} \
    --range ~{variant_range} \
    --minMAC ~{min_mac} \
    --gz \
    --apply-rint \
    --verbose \
    --out ~{outname} 

    >>>
    
    output {

        Array[File]+ regenie_output = glob("~{outname}*.regenie.gz")
        Array[File]+ regenie_log = glob("~{outname}*.log")

    }

    runtime {

		dx_instance_type: "mem2_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/regenie_3.2.1.tar.gz"
    
    }


}

