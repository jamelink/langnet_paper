version 1.1

workflow exome_genotype_variant_filtering {

	input {

		File vcf_file
		File vcf_index
		File subject_list
        File target_regions
		String project_name="my_project"
		Int DP_snv=7
		Int DP_indel=10
		Int GQ=20
		Int average_GQ=35
		Float genotype_missingness=0.10
		Int min_MAC=1
		Float min_AB_indel=0.2
		Float min_AB_snv=0.15

	}

    call extract_subjects_and_target_regions { input: vcf_file = vcf_file, vcf_index = vcf_index, target_regions=target_regions, subject_list = subject_list, project_name = project_name }

    call genotype_filter { input: vcf_name = extract_subjects_and_target_regions.vcf_name, vcf_file_target_regions_extracted_subjects = extract_subjects_and_target_regions.vcf_file_target_regions_extracted_subjects, vcf_index_target_regions_extracted_subjects = extract_subjects_and_target_regions.vcf_index_target_regions_extracted_subjects, project_name = project_name, DP_snv = DP_snv, DP_indel = DP_indel, GQ = GQ }

    call variant_filter_average_gq { input: vcf_name = extract_subjects_and_target_regions.vcf_name, genotype_filtered_vcf = genotype_filter.genotype_filtered_vcf, genotype_filtered_vcf_index = genotype_filter.genotype_filtered_vcf_index, project_name = project_name, average_GQ = average_GQ }

    call variant_filter_genotype_missingness { input: vcf_name = extract_subjects_and_target_regions.vcf_name, genotype_filtered_vcf = genotype_filter.genotype_filtered_vcf, genotype_filtered_vcf_index = genotype_filter.genotype_filtered_vcf_index, project_name = project_name, genotype_missingness = genotype_missingness }

    call variant_filter_allele_count { input: vcf_name = extract_subjects_and_target_regions.vcf_name, genotype_filtered_vcf = genotype_filter.genotype_filtered_vcf, genotype_filtered_vcf_index = genotype_filter.genotype_filtered_vcf_index, project_name = project_name, min_MAC = min_MAC }

    call variant_filter_allele_balance { input: vcf_name = extract_subjects_and_target_regions.vcf_name, genotype_filtered_vcf = genotype_filter.genotype_filtered_vcf, genotype_filtered_vcf_index = genotype_filter.genotype_filtered_vcf_index, variants_fail_mac =  variant_filter_allele_count.variants_fail_mac, project_name = project_name, min_AB_indel = min_AB_indel, min_AB_snv = min_AB_snv }

    call remove_flagged_variants { input: vcf_name = extract_subjects_and_target_regions.vcf_name, genotype_filtered_vcf = genotype_filter.genotype_filtered_vcf, genotype_filtered_vcf_index = genotype_filter.genotype_filtered_vcf_index, variants_fail_average_gq = variant_filter_average_gq.variants_fail_average_gq, variants_fail_genotype_missingness = variant_filter_genotype_missingness.variants_fail_genotype_missingness, variants_fail_mac = variant_filter_allele_count.variants_fail_mac, variants_fail_ab = variant_filter_allele_balance.variants_fail_ab, project_name = project_name }

	output {

		File array_pre_filtering_summary = extract_subjects_and_target_regions.pre_filtering_summary
		File array_variants_fail_average_gq = variant_filter_average_gq.variants_fail_average_gq
		File array_variants_fail_genotype_missingness = variant_filter_genotype_missingness.variants_fail_genotype_missingness
		File array_variants_fail_mac = variant_filter_allele_count.variants_fail_mac
		File array_variants_ref_minor = variant_filter_allele_count.variants_ref_minor
		File array_variants_fail_ab = variant_filter_allele_balance.variants_fail_ab
		File array_variant_filtered_vcf = remove_flagged_variants.variant_filtered_vcf
		File array_variant_filtered_vcf_index = remove_flagged_variants.variant_filtered_vcf_index
		File array_variant_filter_summary = remove_flagged_variants.variant_filter_summary
		File array_variant_filter_stats = remove_flagged_variants.variant_filter_stats

	}

	meta {

        title: "Exome Genotype and Variant Filtering"
        summary: "This workflow runs genotype and variant-level filtering on exome pVCF files."
        description: "This workflow takes as input a pVCF file, associated pVCF index file, subject list, target regions list, and various filtering-specific parameters. It then extracts variants in target regions, individuals in the subject list, runs genotype- and variant-level filtering, and returns filtering summaries and filtered pVCF files with associated index files."
		details: [{
			contactEmail : "Dick.Schijven@mpi.nl",
			contactOrg : "org-psyl"
			}]
    }

    parameter_meta {

        vcf_file: {
            help: "Input pVCF file.",
            patterns: ["*.vcf.gz"]
        }
        vcf_index: {
            help: "Input index file associated with input pVCF file.",
            patterns: ["*.vcf.gz.tbi"]
        }
        subject_list: {
            help: "List of subjects to filter from full data."
        }
		target_regions: {
			help: "List of target genomic regions, with chr, start and end columns."
		}
		project_name: {
			help: "Name of project that is added to output file names."
		}
		DP_snv: {
			help: "Minimum genotype read depth at SNV variant sites."
		}
		DP_indel: {
			help: "Minimum genotype read depth at INDEL variant sites."
		}
		GQ: {
			help: "Minimum genotype quality."
		}
		average_GQ: {
			help: "Minimum variant-level average genotype quality."
		}
		genotype_missingness: {
			help: "Maximum varian-level genotype missingness."
		}
		min_MAC: {
			help: "Minimum variant-level minor allele count."
		}
		min_AB_indel: {
			help: "Minimum variant-level allele balance at INDEL variant sites."
		}
		min_AB_snv: {
			help: "Minimum variant-level allele balance at SNV variant sites."
		}
		
    }


}

task extract_subjects_and_target_regions {

	input {

		File vcf_file
		File vcf_index
		File subject_list
        File target_regions
		String project_name
	
	}

	String vcf_basename = basename( "${vcf_file}", ".vcf.gz")

	command <<<

	set -x -e -o pipefail

	# Extract subjects, extract variant sites in target regions, remove monoallelic variants
	vcftools \
	--gzvcf ~{vcf_file} \
	--keep ~{subject_list} \
    --bed ~{target_regions} \
	--remove-filtered MONOALLELIC \
	--recode -c | bgzip -c > ~{project_name}_~{vcf_basename}_target_regions_extracted_subjects.vcf.gz

	# Make tabix index file
	tabix -p vcf ~{project_name}_~{vcf_basename}_target_regions_extracted_subjects.vcf.gz

	# Report the number of variant sites removed/remaining in/after pre-filtering steps
	printf "1\tPRE_FILTERING_SUMMARY:\n" > ~{project_name}_~{vcf_basename}_pre_filtering_summary.txt
	printf "2\tVCF_file_base_name:\t"~{vcf_basename}"\n" >> ~{project_name}_~{vcf_basename}_pre_filtering_summary.txt
	printf "3\tNumber_of_variant_sites_remaining_before_variant-level_filtering:\t"$( bcftools query -f "%ID\n" ~{project_name}_~{vcf_basename}_target_regions_extracted_subjects.vcf.gz | wc -l | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_basename}_pre_filtering_summary.txt
	printf "\n" >> ~{project_name}_~{vcf_basename}_pre_filtering_summary.txt

	# Remove redundant files
	rm ~{vcf_file}
	rm ~{vcf_index}

	>>>
	
	output {
	
	File vcf_file_target_regions_extracted_subjects = "~{project_name}_~{vcf_basename}_target_regions_extracted_subjects.vcf.gz"
	File vcf_index_target_regions_extracted_subjects = "~{project_name}_~{vcf_basename}_target_regions_extracted_subjects.vcf.gz.tbi"
	String vcf_name = "~{vcf_basename}"
	File pre_filtering_summary = "~{project_name}_~{vcf_basename}_pre_filtering_summary.txt"

	}

	runtime {

		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
		dx_instance_type: "mem3_ssd1_v2_x4"
	
	}

}

task genotype_filter {

	input {

		File vcf_file_target_regions_extracted_subjects
		File vcf_index_target_regions_extracted_subjects
		String vcf_name
		String project_name
		Int? DP_snv
		Int? DP_indel
		Int? GQ

	}

	command <<<
	
	set -x -e -o pipefail

    # Use vcftools to remove monoallelic sites, extract subjects and run min. DP and GQ filters for SNV sites.
    vcftools --gzvcf ~{vcf_file_target_regions_extracted_subjects} \
    --keep-only-indels \
    --minDP ~{DP_indel} \
    --minGQ ~{GQ} \
    --recode -c | bgzip -c > ~{project_name}_~{vcf_name}_genotype_filter_INDEL.vcf.gz

    # Use vcftools to remove monoallelic sites, extract subjects and run min. DP and GQ filters for SNV sites.
    vcftools \
    --gzvcf ~{vcf_file_target_regions_extracted_subjects} \
    --remove-indels \
    --minDP ~{DP_snv} \
    --minGQ ~{GQ} \
    --recode -c | bgzip -c > ~{project_name}_~{vcf_name}_genotype_filter_SNV.vcf.gz

    # Generate a tabix index file for the filtered INDEL and SNV VCFs
    tabix -p vcf ~{project_name}_~{vcf_name}_genotype_filter_INDEL.vcf.gz
    tabix -p vcf ~{project_name}_~{vcf_name}_genotype_filter_SNV.vcf.gz

	# Concatenate the two filtered VCFs
	bcftools concat -a -O z -o ~{project_name}_~{vcf_name}_genotype_filter.vcf.gz ~{project_name}_~{vcf_name}_genotype_filter_INDEL.vcf.gz ~{project_name}_~{vcf_name}_genotype_filter_SNV.vcf.gz

    # Generate a tabix index file for the concatenated filtered VCF file
    tabix -p vcf ~{project_name}_~{vcf_name}_genotype_filter.vcf.gz

	# Remove redundant files
	rm ~{project_name}_~{vcf_name}_genotype_filter_INDEL.vcf.gz
	rm ~{project_name}_~{vcf_name}_genotype_filter_INDEL.vcf.gz.tbi
	rm ~{project_name}_~{vcf_name}_genotype_filter_SNV.vcf.gz
	rm ~{project_name}_~{vcf_name}_genotype_filter_SNV.vcf.gz.tbi

	>>>

	output {

		File genotype_filtered_vcf = "~{project_name}_~{vcf_name}_genotype_filter.vcf.gz"
		File genotype_filtered_vcf_index = "~{project_name}_~{vcf_name}_genotype_filter.vcf.gz.tbi"

	}

	runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"

	}

}

task variant_filter_average_gq {

	input {

		File genotype_filtered_vcf
		File genotype_filtered_vcf_index
		String vcf_name
		String project_name
		Int? average_GQ

	}

	command <<<

	set -x -e -o pipefail

    # Query all the individual genotype quality (GQ) values from the VCF file
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GQ\t]\n' ~{genotype_filtered_vcf} > ~{vcf_name}_stats_GQ.txt

    # Compute the average genotype quality for each variant site
    awk '{sum = 0; for (i = 6; i <= NF; i++) sum += $i; sum /= (NF-5); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"sum}' ~{vcf_name}_stats_GQ.txt > ~{vcf_name}_stats_average_GQ.txt

    # Get variant sites for which the average GQ is below the set threshold
    awk -v average_gq=~{average_GQ} '$6 < average_gq {print $3}' ~{vcf_name}_stats_average_GQ.txt > ~{project_name}_~{vcf_name}_variants_fail_average_gq.txt

    # If no variant sites are below the average GQ threshold, make an empty file
    if [ ! -f ~{project_name}_~{vcf_name}_variants_fail_average_gq.txt ]; then
        echo "No variants flagged for failing average GQ."
        touch ~{project_name}_~{vcf_name}_variants_fail_average_gq.txt
    fi

	# Remove redundant files
	rm ~{vcf_name}_stats_GQ.txt
	rm ~{vcf_name}_stats_average_GQ.txt

	>>>

	output {

		File variants_fail_average_gq = "~{project_name}_~{vcf_name}_variants_fail_average_gq.txt"
	
	}

	runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
	
	}

}

task variant_filter_genotype_missingness {

	input {

		File genotype_filtered_vcf
		File genotype_filtered_vcf_index
		String vcf_name
		String project_name
		Float? genotype_missingness

	}

	command <<<
	
	set -x -e -o pipefail

    # Make a missingness report for variant sites
    vcftools \
    --gzvcf ~{genotype_filtered_vcf} \
    --missing-site \
    --out ~{vcf_name}_missingness_report

    # Get variant for which genotype missingness is above the threshold
    awk -v genotype_missing=~{genotype_missingness} '$6 >= genotype_missing {print $1"\t"$2}' ~{vcf_name}_missingness_report.lmiss | sed '1d' > ~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt
    
	# If no variant sites are above the genotype missingness threshold, make an empty file
	# else
	# Query the genotype filtered vcf file for the IDs of the variants that fail missingness
    if [ ! -f ~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt ]; then
        echo "No variants flagged for failing genotype missingness."
        touch ~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt
    else	
	bcftools query -T ~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt -f '%ID\n' ~{genotype_filtered_vcf} > tmp && mv tmp ~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt
	fi

	# Remove redundant files
	rm ~{vcf_name}_missingness_report.lmiss

	>>>

	output {

		File variants_fail_genotype_missingness = "~{project_name}_~{vcf_name}_variants_fail_genotype_missingness.txt"
	
	}

	runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
	
	}

}

task variant_filter_allele_count {

	input {

		File genotype_filtered_vcf
		File genotype_filtered_vcf_index
		String vcf_name
		String project_name
		Int? min_MAC

	}

	command <<<
	
	set -x -e -o pipefail

	# Make an allele counts report for variant sites
	vcftools \
	--gzvcf ~{genotype_filtered_vcf} \
	--counts2 \
	--out  ~{vcf_name}_mac

	# Report variants where the reference allele is minor
	awk -v OFS="\t" '$5 < $4/2 {print $1,$2,$3,$4,$5,$4-$5}' ~{vcf_name}_mac.frq.count > ~{project_name}_~{vcf_name}_variants_ref_is_minor.txt
	printf "CHR\tPOS\tN_ALLELE\tN_REF\tN_ALT\n" > header.txt
	cat header.txt ~{project_name}_~{vcf_name}_variants_ref_is_minor.txt > tmp && mv tmp ~{project_name}_~{vcf_name}_variants_ref_is_minor.txt

	# Remove redundant files
	rm header.txt

	# Remove header from the counts report, replace all : by \t
	sed '1d' ~{vcf_name}_mac.frq.count > ~{vcf_name}_mac.frq.count.tmp

	# Run through the lines of the counts report
	while read line; do

	# Save line as a temporary file
	echo ${line} > ~{vcf_name}_mac.frq.count.tmp.line

	# Starting at column 5, print the MAC value in each column.
	# Count the number of lines where the MAC is lower than the set threshold.
	# If that number equals the total number of alleles on that line - 1 (because the ref allele is always expected to have MAC > 0), flag the variant for removal.
	# For multi-allelic variants, these are only flagged when ALL alternate alleles have a MAC below the threshold.
	if [[ $( awk '{for (i=5; i<=NF; i++) {print $i}}'  ~{vcf_name}_mac.frq.count.tmp.line | awk -v mac=~{min_MAC} '$1 < mac' | wc -l ) -eq $(( $( awk '{for (i=5; i<=NF; i++) {print $i}}'  ~{vcf_name}_mac.frq.count.tmp.line | wc -l )-1 )) ]]; then

	awk '{print $1"\t"$2}'  ~{vcf_name}_mac.frq.count.tmp.line >> ~{project_name}_~{vcf_name}_variants_fail_mac.txt

	fi

	# Starting at column 5, print the MAC value in each column.
	# Count the number of lines where the MAC is equal to 0.
	# If that number equals the total number of alleles on that line, this means that all genotypes are missing, thus MAC = 0, and flag the variant for removal.
	if [[ $( awk '{for (i=5; i<=NF; i++) {print $i}}'  ~{vcf_name}_mac.frq.count.tmp.line | awk -v mac=~{min_MAC} '$1 == 0' | wc -l ) -eq $( awk '{for (i=5; i<=NF; i++) {print $i}}'  ~{vcf_name}_mac.frq.count.tmp.line | wc -l ) ]]; then

	awk '{print $1"\t"$2}'  ~{vcf_name}_mac.frq.count.tmp.line >> ~{project_name}_~{vcf_name}_variants_fail_mac.txt

	fi

	# Remove the line file
	rm ~{vcf_name}_mac.frq.count.tmp.line

	done < ~{vcf_name}_mac.frq.count.tmp

	# If no variant sites are below the MAC threshold, make an empty file
		if [ ! -f ~{project_name}_~{vcf_name}_variants_fail_mac.txt ]; then
	echo "No variants flagged for failing MAC."
	touch ~{project_name}_~{vcf_name}_variants_fail_mac.txt
	else
	# else
	# Query the genotype filtered vcf file for the IDs of the variants that fail missingness
	# First remove variants that do not pass filter from the report of minor reference alleles
	awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' ~{project_name}_~{vcf_name}_variants_fail_mac.txt ~{project_name}_~{vcf_name}_variants_ref_is_minor.txt > tmp && mv tmp ~{project_name}_~{vcf_name}_variants_ref_is_minor.txt
	bcftools query -T ~{project_name}_~{vcf_name}_variants_fail_mac.txt -f '%ID\n' ~{genotype_filtered_vcf} > tmp && mv tmp ~{project_name}_~{vcf_name}_variants_fail_mac.txt
	fi

	

	# Remove redundant files
	rm ~{vcf_name}_mac.frq.count
	rm ~{vcf_name}_mac.frq.count.tmp

	>>>

	output {

		File variants_fail_mac = "~{project_name}_~{vcf_name}_variants_fail_mac.txt"
		File variants_ref_minor = "~{project_name}_~{vcf_name}_variants_ref_is_minor.txt"
	
	}

	runtime {

		dx_instance_type: "mem1_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
	
	}

}

task variant_filter_allele_balance {

	input {

		File genotype_filtered_vcf
		File genotype_filtered_vcf_index
		File variants_fail_mac
		String vcf_name
		String project_name
		Float? min_AB_indel
		Float? min_AB_snv

	}

	command <<<
	
	set -x -e -o pipefail

	# Query variant sites including heterozygous genotypes with alternate allele, and excluding homozygous genotypes with alternate allele
	bcftools view -i 'GT="het"' ~{genotype_filtered_vcf} | bcftools query -e 'GT="AA"' -f '%ID\n' > ~{vcf_name}_ABcheck.txt

	# Extract INDEL sites for AB check
	vcftools \
	--gzvcf ~{genotype_filtered_vcf} \
	--snps ~{vcf_name}_ABcheck.txt \
	--exclude ~{variants_fail_mac} \
	--keep-only-indels \
	--recode -c | bgzip -c > ~{vcf_name}_ABcheck_INDEL.vcf.gz

	# Generate tabix index file for AB check INDEL VCF
	tabix -p vcf ~{vcf_name}_ABcheck_INDEL.vcf.gz

	# Extract SNV sites for AB check
	vcftools \
	--gzvcf ~{genotype_filtered_vcf} \
	--snps ~{vcf_name}_ABcheck.txt \
	--exclude ~{variants_fail_mac} \
	--remove-indels \
	--recode -c | bgzip -c > ~{vcf_name}_ABcheck_SNV.vcf.gz

	# Generate tabix index file for AB check SNV VCF
	tabix -p vcf ~{vcf_name}_ABcheck_SNV.vcf.gz

	# Run AB check for INDEL VCF
	gatk --java-options "-Xmx12G" FilterVcf \
	-I ~{vcf_name}_ABcheck_INDEL.vcf.gz \
	-O ~{vcf_name}_ABcheck_INDEL_checked.vcf.gz \
	--MIN_AB ~{min_AB_indel}

	# Get variant sites that do not pass AB check in INDEL VCF
	bcftools query -e 'FILTER="PASS"' -f '%ID\n' ~{vcf_name}_ABcheck_INDEL_checked.vcf.gz > ~{vcf_name}_remove_indel_sites_AB.txt

	# Run AB check for SNV VCF
	gatk --java-options "-Xmx12G" FilterVcf \
	-I ~{vcf_name}_ABcheck_SNV.vcf.gz \
	-O ~{vcf_name}_ABcheck_SNV_checked.vcf.gz \
	--MIN_AB ~{min_AB_snv}

	# Get variant sites that do not pass AB check in SNV VCF
	bcftools query -e 'FILTER="PASS"' -f '%ID\n' ~{vcf_name}_ABcheck_SNV_checked.vcf.gz > ~{vcf_name}_remove_snv_sites_AB.txt

	# Concatenate variant sites that fail AB filter into one file
	cat ~{vcf_name}_remove_indel_sites_AB.txt ~{vcf_name}_remove_snv_sites_AB.txt | sort | uniq > ~{project_name}_~{vcf_name}_variants_fail_ab.txt

	# Remove redundant files
	rm ~{vcf_name}_ABcheck.txt
	rm ~{vcf_name}_ABcheck_INDEL.vcf.gz
	rm ~{vcf_name}_ABcheck_INDEL.vcf.gz.tbi
	rm ~{vcf_name}_ABcheck_SNV.vcf.gz
	rm ~{vcf_name}_ABcheck_SNV.vcf.gz.tbi
	rm ~{vcf_name}_remove_indel_sites_AB.txt
	rm ~{vcf_name}_remove_snv_sites_AB.txt

	>>>
	
	output {

		File variants_fail_ab = "~{project_name}_~{vcf_name}_variants_fail_ab.txt"
	
	}

	runtime {

		dx_instance_type: "mem2_ssd1_v2_x4"
		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
	
	}

}

task remove_flagged_variants {
	
	input {

		File genotype_filtered_vcf
		File genotype_filtered_vcf_index
		String vcf_name
		File variants_fail_average_gq
		File variants_fail_genotype_missingness
		File variants_fail_mac
		File variants_fail_ab
		String project_name
	
	}
	
	command <<<
	
	set -x -e -o pipefail

	# Concatenate and summarize the variants flagged for removal
	cat ~{variants_fail_average_gq} ~{variants_fail_mac} ~{variants_fail_genotype_missingness} ~{variants_fail_ab} | sort | uniq > ~{vcf_name}_variants_fail_filtering_all.txt

	# Report the numbers of sites flagged for removal based on filtering criteria
	printf "1\tVARIANT_FILTERING_SUMMARY:\n" > ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "2\tVCF_file_base_name:\t"~{vcf_name}"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "3\tAverage_genotype_quality_across_genotypes:\t"$( wc -l ~{variants_fail_average_gq} | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "4\tMissing_genotype_frequency:\t"$( wc -l ~{variants_fail_genotype_missingness} | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "5\tMinor_Allele_Count_(MAC):\t"$( wc -l ~{variants_fail_mac} | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "6\tAllele_Balance_(AB)_after_removing_sites_with_low_MAC:\t"$( wc -l ~{variants_fail_ab} | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "7\tTotal_number_of_unique_sites_flagged_for_removal:\t"$( wc -l ~{vcf_name}_variants_fail_filtering_all.txt | awk '{print $1}' )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt

	# Remove flagged sites from the data
	vcftools \
	--gzvcf ~{genotype_filtered_vcf} \
	--exclude ~{vcf_name}_variants_fail_filtering_all.txt \
	--recode -c | bgzip -c > ~{project_name}_~{vcf_name}_variant_filter.vcf.gz

	tabix -p vcf ~{project_name}_~{vcf_name}_variant_filter.vcf.gz

	# Calculate TsTv ratio pre-filtering
	vcftools \
	--gzvcf ~{genotype_filtered_vcf} \
	--TsTv-summary \
	--out ~{vcf_name}_pre_filtering

	# Report the TsTv ratio pre-filtering in the filtering summary file
	printf "8\tTs_pre-filtering:\t"$( grep Ts ~{vcf_name}_pre_filtering.TsTv.summary | awk '{print $2}')"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "9\tTv_pre-filtering:\t"$( grep Tv ~{vcf_name}_pre_filtering.TsTv.summary | awk '{print $2}')"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt

	# Calculate TsTv ratio post-filtering
	vcftools \
	--gzvcf ~{project_name}_~{vcf_name}_variant_filter.vcf.gz \
	--TsTv-summary \
	--out ~{vcf_name}_post_filtering

	# Report the TsTv ratio post-filtering in the filtering summary file
	printf "10\tTs_post-filtering:\t"$( grep Ts ~{vcf_name}_post_filtering.TsTv.summary | awk '{print $2}')"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt
	printf "11\tTv_post-filtering:\t"$( grep Tv ~{vcf_name}_post_filtering.TsTv.summary | awk '{print $2}')"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt

	# Output statistics of the VCF file
	bcftools stats -d 0,40,1 ~{project_name}_~{vcf_name}_variant_filter.vcf.gz > ~{project_name}_~{vcf_name}_variant_filter_stats.txt

	printf "12\tTotal_number_of_variant_sites_remaining:\t"$( awk -F"\t" '$1 == "SN" && $3 == "number of records:" {print $4}' ~{project_name}_~{vcf_name}_variant_filter_stats.txt )"\n" >> ~{project_name}_~{vcf_name}_variant_filter_summary.txt

	>>>

	output {

		File variant_filtered_vcf = "~{project_name}_~{vcf_name}_variant_filter.vcf.gz"
    	File variant_filtered_vcf_index = "~{project_name}_~{vcf_name}_variant_filter.vcf.gz.tbi"
    	File variant_filter_summary = "~{project_name}_~{vcf_name}_variant_filter_summary.txt"
		File variant_filter_stats = "~{project_name}_~{vcf_name}_variant_filter_stats.txt"

	}

	runtime {

		docker: "dx://AICHA_rs_conn:/workflows/docker_images/exome_filtering_1.0.tar.gz"
		dx_instance_type: "mem1_ssd1_v2_x4"

	}


}

