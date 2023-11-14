version 1.1

workflow exome_merge_plink_files {

	input {

		Array[File]+ plink_bed_files
		Array[File]+ plink_bim_files
        Array[File]+ plink_fam_files
        String out_prefix
        Boolean split_x = false
        File? individuals_sex

	}


    call merge_plink_files { input: plink_bed_files = plink_bed_files, plink_bim_files = plink_bim_files, plink_fam_files = plink_fam_files, out_prefix = out_prefix, split_x = split_x, individuals_sex = individuals_sex }


    call convert_to_bgen { input: merged_plink_bed_file = merge_plink_files.merged_plink_bed_file, merged_plink_bim_file = merge_plink_files.merged_plink_bim_file, merged_plink_fam_file = merge_plink_files.merged_plink_fam_file, merged_plink_bed_file_PAR1 = merge_plink_files.merged_plink_bed_file_PAR1, merged_plink_bim_file_PAR1 = merge_plink_files.merged_plink_bim_file_PAR1, merged_plink_fam_file_PAR1 = merge_plink_files.merged_plink_fam_file_PAR1, merged_plink_bed_file_PAR2 = merge_plink_files.merged_plink_bed_file_PAR2, merged_plink_bim_file_PAR2 = merge_plink_files.merged_plink_bim_file_PAR2, merged_plink_fam_file_PAR2 = merge_plink_files.merged_plink_fam_file_PAR2, merged_plink_bed_file_XnoPAR = merge_plink_files.merged_plink_bed_file_XnoPAR, merged_plink_bim_file_XnoPAR = merge_plink_files.merged_plink_bim_file_XnoPAR, merged_plink_fam_file_XnoPAR = merge_plink_files.merged_plink_fam_file_XnoPAR, merged_plink_bed_file_Xclean = merge_plink_files.merged_plink_bed_file_Xclean, merged_plink_bim_file_Xclean = merge_plink_files.merged_plink_bim_file_Xclean, merged_plink_fam_file_Xclean = merge_plink_files.merged_plink_fam_file_Xclean, split_x = split_x, out_prefix = out_prefix }


    meta {

        title: "Exome merge plink files"
        summary: "This workflow merges exome plink files."
        description: "This workflow takes a input an array of plink .bed, .bim and .fam files and merges these binary filesets. Optionally it outputs a merged bgen fileset."
    
    }

    parameter_meta {

        plink_bed_files: {
            help: "An array of input pVCF files.",
            patterns: ["*.bed"]
        }
        plink_bim_files: {
            help: "An array of input pVCF index files. The order of these should correspond with the associated input pVCF files.",
            patterns: ["*.bim"]
        }
        plink_fam_files: {
            help: "An array of input pVCF index files. The order of these should correspond with the associated input pVCF files.",
            patterns: ["*.fam"]
        }
        out_prefix: {
            help: "Prefix of output files."
        }
        split_x: {
            help: "Whether or not to split X chromosome data into PAR1, PAR2 and XnoPAR (the X chromosome where het. haploid genotypes for males are set to missing)."
        }
    
    }

}

task merge_plink_files {

    input {

        Array[File]+ plink_bed_files
        Array[File]+ plink_bim_files
        Array[File]+ plink_fam_files
        String out_prefix
        Boolean split_x
        File? individuals_sex

    }

    Int plink_files = length(plink_bed_files)

    command <<<

    set -x -e -o pipefail

    # Move input files to a new directory
    mkdir /input_files

    mv ~{sep(" ", plink_bed_files)} /input_files/
    mv ~{sep(" ", plink_bim_files)} /input_files/
    mv ~{sep(" ", plink_fam_files)} /input_files/

    # Make a list of files that need to be merged
    ls /input_files/*.bed > merge_bed_files.txt
    ls /input_files/*.bim > merge_bim_files.txt
    ls /input_files/*.fam > merge_fam_files.txt

    first_file=$( head -1 merge_bed_files.txt | awk -F".bed" '{print $1}' )

    paste merge_bed_files.txt merge_bim_files.txt merge_fam_files.txt | sed '1d' > merge_list.txt

    # Merge the first file with the files in the input list
    plink \
    --bed ${first_file}.bed \
    --bim ${first_file}.bim \
    --fam ${first_file}.fam \
    --merge-list merge_list.txt \
    --make-bed \
    --out ~{out_prefix}

    if ~{split_x}; then

        # Split off the PAR1 region, leave this as is
        plink \
        --bfile ~{out_prefix} \
        --chr 23 \
        --update-sex ~{individuals_sex} \
        --to-bp 2781479 \
        --make-bed \
        --out ~{out_prefix}_PAR1

        # Split off the PAR2 region, leave this as is
        plink \
        --bfile ~{out_prefix} \
        --chr 23 \
        --update-sex ~{individuals_sex} \
        --from-bp 155701383 \
        --make-bed \
        --out ~{out_prefix}_PAR2

        # Split off the rest of chromosome X, set heterozygous haploid genotypes to missing
        plink \
        --chr 23 \
        --bfile ~{out_prefix} \
        --update-sex ~{individuals_sex} \
        --from-bp 2781480 \
        --to-bp 155701382 \
        --set-hh-missing \
        --make-bed \
        --out ~{out_prefix}_XnoPAR

        # Make a missingness report for hh genotypes, which allows to exlucde variants later based on a threshold
        plink \
        --bfile ~{out_prefix}_XnoPAR \
        --extract ~{out_prefix}_XnoPAR.hh \
        --missing \
        --out ~{out_prefix}_XnoPAR

        # Merge the PAR1, PAR2 and XnoPAR regions
        if [ -f ~{out_prefix}_PAR2.bed ]; then
        
            printf "~{out_prefix}_XnoPAR.bed\t~{out_prefix}_XnoPAR.bim\t~{out_prefix}_XnoPAR.fam\n~{out_prefix}_PAR2.bed\t~{out_prefix}_PAR2.bim\t~{out_prefix}_PAR2.fam\n" > merge_list_cX.txt
        
        else

            printf "~{out_prefix}_XnoPAR.bed\t~{out_prefix}_XnoPAR.bim\t~{out_prefix}_XnoPAR.fam\n" > merge_list_cX.txt

        fi

        if [ -f ~{out_prefix}_PAR1.bed ]; then
        
        plink \
        --bfile ~{out_prefix}_PAR1 \
        --merge-list merge_list_cX.txt \
        --make-bed \
        --out ~{out_prefix}_Xclean

        else

        mv ~{out_prefix}_XnoPAR.bed ~{out_prefix}_Xclean.bed
        mv ~{out_prefix}_XnoPAR.bim ~{out_prefix}_Xclean.bim
        mv ~{out_prefix}_XnoPAR.fam ~{out_prefix}_Xclean.fam

        fi

    fi

    >>>
    
    output {

        File merged_plink_bed_file = "~{out_prefix}.bed"
        File merged_plink_bim_file = "~{out_prefix}.bim"
        File merged_plink_fam_file = "~{out_prefix}.fam"
        File? merged_plink_bed_file_PAR1 = "~{out_prefix}_PAR1.bed"
        File? merged_plink_bim_file_PAR1 = "~{out_prefix}_PAR1.bim"
        File? merged_plink_fam_file_PAR1 = "~{out_prefix}_PAR1.fam"
        File? merged_plink_bed_file_PAR2 = "~{out_prefix}_PAR2.bed"
        File? merged_plink_bim_file_PAR2 = "~{out_prefix}_PAR2.bim"
        File? merged_plink_fam_file_PAR2 = "~{out_prefix}_PAR2.fam"
        File? merged_plink_bed_file_XnoPAR = "~{out_prefix}_XnoPAR.bed"
        File? merged_plink_bim_file_XnoPAR = "~{out_prefix}_XnoPAR.bim"
        File? merged_plink_fam_file_XnoPAR = "~{out_prefix}_XnoPAR.fam"
        File? chrX_missing_report = "~{out_prefix}_XnoPAR.lmiss"
        File? merged_plink_bed_file_Xclean = "~{out_prefix}_Xclean.bed"
        File? merged_plink_bim_file_Xclean = "~{out_prefix}_Xclean.bim"
        File? merged_plink_fam_file_Xclean = "~{out_prefix}_Xclean.fam"

    }

    runtime {

		dx_instance_type: "mem1_ssd2_v2_x8"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"

	}


}

task convert_to_bgen {

    input {
        
        File merged_plink_bed_file
        File merged_plink_bim_file
        File merged_plink_fam_file
        File? merged_plink_bed_file_PAR1
        File? merged_plink_bim_file_PAR1
        File? merged_plink_fam_file_PAR1
        File? merged_plink_bed_file_PAR2
        File? merged_plink_bim_file_PAR2
        File? merged_plink_fam_file_PAR2
        File? merged_plink_bed_file_XnoPAR
        File? merged_plink_bim_file_XnoPAR
        File? merged_plink_fam_file_XnoPAR
        File? merged_plink_bed_file_Xclean
        File? merged_plink_bim_file_Xclean
        File? merged_plink_fam_file_Xclean
        String out_prefix
        Boolean split_x
       
    }

    command <<<

    set -x -e -o pipefail

    # Convert the plink binary fileset to bgen 1.2 8-bits format
    plink2 \
    --bed ~{merged_plink_bed_file} \
    --bim ~{merged_plink_bim_file} \
    --fam ~{merged_plink_fam_file} \
    --export bgen-1.2 bits=8 ref-first \
    --out ~{out_prefix}

    # Make bgen index file
    bgenix -g ~{out_prefix}.bgen \
    -index \
    -clobber


   if ~{split_x}; then

    # Convert the plink binary fileset to bgen 1.2 8-bits format
    plink2 \
    --bed ~{merged_plink_bed_file_PAR1} \
    --bim ~{merged_plink_bim_file_PAR1} \
    --fam ~{merged_plink_fam_file_PAR1} \
    --export bgen-1.2 bits=8 ref-first \
    --out ~{out_prefix}_PAR1

    # Make bgen index file
    bgenix -g ~{out_prefix}_PAR1.bgen \
    -index \
    -clobber

    # Convert the plink binary fileset to bgen 1.2 8-bits format
    plink2 \
    --bed ~{merged_plink_bed_file_PAR2} \
    --bim ~{merged_plink_bim_file_PAR2} \
    --fam ~{merged_plink_fam_file_PAR2} \
    --export bgen-1.2 bits=8 ref-first \
    --out ~{out_prefix}_PAR2

    # Make bgen index file
    bgenix -g ~{out_prefix}_PAR2.bgen \
    -index \
    -clobber

    # Convert the plink binary fileset to bgen 1.2 8-bits format
    plink2 \
    --bed ~{merged_plink_bed_file_XnoPAR} \
    --bim ~{merged_plink_bim_file_XnoPAR} \
    --fam ~{merged_plink_fam_file_XnoPAR} \
    --export bgen-1.2 bits=8 ref-first \
    --out ~{out_prefix}_XnoPAR

    # Make bgen index file
    bgenix -g ~{out_prefix}_XnoPAR.bgen \
    -index \
    -clobber

    # Convert the plink binary fileset to bgen 1.2 8-bits format
    plink2 \
    --bed ~{merged_plink_bed_file_Xclean} \
    --bim ~{merged_plink_bim_file_Xclean} \
    --fam ~{merged_plink_fam_file_Xclean} \
    --export bgen-1.2 bits=8 ref-first \
    --out ~{out_prefix}_Xclean

    # Make bgen index file
    bgenix -g ~{out_prefix}_Xclean.bgen \
    -index \
    -clobber

    fi

    >>>

    output {

        File bgen_file = "~{out_prefix}.bgen"
        File bgen_sample = "~{out_prefix}.sample"
        File bgen_index = "~{out_prefix}.bgen.bgi"
        File? bgen_file_PAR1 = "~{out_prefix}_PAR1.bgen"
        File? bgen_sample_PAR1 = "~{out_prefix}_PAR1.sample"
        File? bgen_index_PAR1 = "~{out_prefix}_PAR1.bgen.bgi"
        File? bgen_file_PAR2 = "~{out_prefix}_PAR2.bgen"
        File? bgen_sample_PAR2 = "~{out_prefix}_PAR2.sample"
        File? bgen_index_PAR2 = "~{out_prefix}_PAR2.bgen.bgi"
        File? bgen_file_XnoPAR = "~{out_prefix}_XnoPAR.bgen"
        File? bgen_sample_XnoPAR = "~{out_prefix}_XnoPAR.sample"
        File? bgen_index_XnoPAR = "~{out_prefix}_XnoPAR.bgen.bgi"
        File? bgen_file_Xclean = "~{out_prefix}_Xclean.bgen"
        File? bgen_sample_Xclean = "~{out_prefix}_Xclean.sample"
        File? bgen_index_Xclean = "~{out_prefix}_Xclean.bgen.bgi"
        
    }

    runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"

	}

}