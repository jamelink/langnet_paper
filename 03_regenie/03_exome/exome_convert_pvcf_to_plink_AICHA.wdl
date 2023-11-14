version 1.1

workflow exome_convert_pvcf_to_plink {

	input {

		File vcf_file
		File vcf_index

	}

    call remove_multiallelic_variants { input: vcf_file = vcf_file, vcf_index = vcf_index }

    call convert_to_plink { input: biallelic_only_vcf_file = remove_multiallelic_variants.biallelic_only_vcf_file, vcf_name = remove_multiallelic_variants.vcf_name }

    output {

    File plink_bed_file = convert_to_plink.plink_bed_file
    File plink_bim_file = convert_to_plink.plink_bim_file
    File plink_fam_file = convert_to_plink.plink_fam_file
    File conversion_summary = remove_multiallelic_variants.conversion_summary

    }

    meta {

        title: "Convert exome pVCF to plink"
        summary: "This workflow converts a pVCF file to plink file format."
        description: "This workflow takes a input a pVCF file and associated pVCF index. It then removes multiallelic variant sites and converts the pVCF file to plink binary format."
    
    }

    parameter_meta {

        vcf_file: {
            help: "An array of input pVCF files.",
            patterns: ["*.vcf.gz"]
        }
        vcf_index: {
            help: "An array of input pVCF index files. The order of these should correspond with the associated input pVCF files.",
            patterns: ["*.vcf.gz.tbi"]
        }
    
    }

}

task remove_multiallelic_variants {

    input {

        File vcf_file
        File vcf_index
    
    }

    String vcf_basename = basename( "${vcf_file}", ".vcf.gz")

    command <<<

    set -x -e -o pipefail
    
    # Print only the biallelic variants in the vcf file
    bcftools view --max-alleles 2 -O z -o ~{vcf_basename}.biallelic.vcf.gz ~{vcf_file}

    # Make a short summary containing the number of removed multiallelic variant sites
    printf "1\tCONVERSION_TO_ANALYSIS_READY_FORMAT:\n" > ~{vcf_basename}_conversion_summary.txt
	printf "2\tVCF_file_base_name:\t"~{vcf_basename}"\n" >> ~{vcf_basename}_conversion_summary.txt
	printf "3\tNumber_of_multiallelic_variant_sites_removed:\t"$( bcftools view --min-alleles 3 ~{vcf_file} | bcftools query -f '%ID\n' | wc -l )"\n" >> ~{vcf_basename}_conversion_summary.txt

    >>>

    output {

        File biallelic_only_vcf_file = "~{vcf_basename}.biallelic.vcf.gz"
        File conversion_summary = "~{vcf_basename}_conversion_summary.txt"
        String vcf_name = "~{vcf_basename}"

    }
    
    runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"

	}

}


task convert_to_plink {

    input {

       File biallelic_only_vcf_file
       String vcf_name
    
    }

    command <<<

    set -x -e -o pipefail

    plink \
    --vcf ~{biallelic_only_vcf_file} \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --double-id \
    --allow-extra-chr 0 \
    --make-bed \
    --vcf-half-call m \
    --out ~{vcf_name}

    >>>

    output {

        File plink_bed_file = "~{vcf_name}.bed"
        File plink_bim_file = "~{vcf_name}.bim"
        File plink_fam_file = "~{vcf_name}.fam"

    }

   	runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"

	}

}