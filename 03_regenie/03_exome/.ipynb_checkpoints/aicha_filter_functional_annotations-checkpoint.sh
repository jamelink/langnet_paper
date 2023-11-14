#!/bin/sh
#$ -N filter_func_annot
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

chr=$1

###############################
## SET DIRECTORIES AND FILES ##
###############################

# RunName: This is the same run name that was specified for the genotype and variant filtering
RunName=AICHA_rs_conn

# Location where plink .bim files of the filtered variant output are stored
PLINKDir=/data/clusterfs/lag/users/jitame/SENT_CORE/geno/regenie/exome/functional_annotation
# Output directory where functional annotation filtering will be stored
OutDir=/data/clusterfs/lag/users/jitame/SENT_CORE/geno/regenie/exome/functional_annotation/selected_variants_v1

# File with canonical transcripts (first column should be Ensembl gene ID, second column should be Ensembl transcript ID)
CanonicalTranscripts=/data/clusterfs/lag/users/jitame/SENT_CORE/geno/regenie/exome/functional_annotation/MANE.GRCh38.v1.0.gene_transcript.txt

# CADD score thresholds for strict and broad filters
CADD_strict=20
CADD_broad=1

# File containing the chromosome-block files (do not change)
VCFblocks=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/helper_files/pvcf_blocks.txt

# Directory where the raw annotation files are stored (do not change)
AnnotDir=/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/exome/exome_release_final/functional_annotation/annotated_variants

# Directory where extracted dbNSFP scores are stored (do not change)
dbnsfpDir=/data/clusterfs/lag/users/jitame/SENT_CORE/geno/regenie/exome/functional_annotation/dbnsfp_cadd

# Get the blocks for this chromosome
blocks=$( awk -F"\t" -v chr=${chr} '$2 == chr {print $3}' ${VCFblocks} )

###############################################
## SELECT VARIANTS WITH HIGH PUTATIVE IMPACT ##
###############################################

echo "SELECT VARIANTS WITH HIGH PUTATIVE IMPACT"

for block in ${blocks}; do

printf ${block}" "

# Select protein-coding variants with a high putative impact, and only variants that are not in the 5% tail ends of the protein, then print column 3 (variant ID), column 10 (ENSEMBL gene ID) and column 12 (ENSEMBL transcript ID)
awk -F "\t" '$8 == "HIGH" && $13 == "protein_coding" && ($21 > $22*0.05) && ($21 < $22*0.95) {print $3"\t"$10"\t"$12}' ${AnnotDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.tab | awk -F"." '{print $1}' >> ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_high.txt

done

cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_high.txt > ${OutDir}/${RunName}_c${chr}_variants_impact_high.txt
rm ${OutDir}/${RunName}_c${chr}_b*_variants_impact_high.txt

# Canonical transcripts
awk 'FNR==NR{seen[$0];next}($2"\t"$3 in seen){print}' ${CanonicalTranscripts} ${OutDir}/${RunName}_c${chr}_variants_impact_high.txt | awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim - | awk '{print $1"\t"$2"\thigh_canonical"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt

# Alternative transcripts
awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_variants_impact_high.txt | awk '{print $1"\t"$2"\thigh_alt"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt

# Exclude variant-gene combinations that are already in the canonical list
awk 'NR==FNR{a[$1"\t"2];next} !($1"\t"2 in a)' ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt

#############################################
## SELECT ADDITIONAL VARIANTS FOR THE SETS ##
#############################################

printf "\n\nSELECT ADDITIONAL VARIANTS FOR THE SETS\n"

for block in ${blocks}; do

    printf ${block}" "

    ### MODERATE IMPACT VARIANTS

    awk -F"\t" '$8 == "MODERATE" && $13 == "protein_coding" {print $3"\t"$10"\t"$12}' ${AnnotDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.tab | awk -F"." '{print $1}' > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt
    
    col_n_cadd_phred=$( awk -v RS='\t' '/dbNSFP_CADD_phred/{print NR; exit}' ${dbnsfpDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.dbnsfp.cadd.singleline.tab )

    awk -F"\t" -v col_n_cadd_phred=${col_n_cadd_phred} -v cs=${CADD_strict} '$col_n_cadd_phred >= cs && $3 != "ID" {print $3"\t"$6"\t"$7}' ${dbnsfpDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.dbnsfp.cadd.singleline.tab > ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_strict.txt
    awk -F"\t" -v col_n_cadd_phred=${col_n_cadd_phred} -v cb=${CADD_broad} -v cs=${CADD_strict} '$col_n_cadd_phred >= cb && $col_n_cadd_phred < cs  && $3 != "ID" {print $3"\t"$6"\t"$7}' ${dbnsfpDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.dbnsfp.cadd.singleline.tab > ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_broad.txt

    # Strict: CADD >= 20 
    # Canonical transcripts ###OK!
    awk 'FNR==NR{seen[$0];next}($2"\t"$3 in seen){print}' ${CanonicalTranscripts} ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt | awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_strict.txt - | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_canonical.txt
    awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_canonical.txt > ${OutDir}/c${chr}_tmp && mv ${OutDir}/c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_canonical.txt
    
    # Alternative transcripts
    awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_strict.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_alt_transcripts.txt
    awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_alt_transcripts.txt
    awk 'NR==FNR{a[$0];next} !($0 in a)' ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_canonical.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_strict_alt_transcripts.txt
    
    # Broad: CADD >= 1 && CADD < 20
    # Canonical transcripts
    awk 'FNR==NR{seen[$0];next}($2"\t"$3 in seen){print}' ${CanonicalTranscripts} ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt | awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_broad.txt - | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_canonical.txt
    awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_canonical.txt > ${OutDir}/c${chr}_tmp && mv ${OutDir}/c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_canonical.txt

    # Alternative transcripts
    awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_broad.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_alt_transcripts.txt
    awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_alt_transcripts.txt
    awk 'NR==FNR{a[$0];next} !($0 in a)' ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_canonical.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate_broad_alt_transcripts.txt
    

    ### MODIFIER IMPACT VARIANTS
    awk -F"\t" '$8 == "MODIFIER" && $13 == "protein_coding" {print $3"\t"$10"\t"$12}' ${AnnotDir}/c${chr}/ukb23157_c${chr}_b${block}_v1_site_only.snpeff.tab | awk -F"." '{print $1}' > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier.txt

    #Broad: CADD >= 1
    awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_strict.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier.txt | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_strict.txt
    awk 'FNR==NR{seen[$0];next}($0 in seen){print}' ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_broad.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier.txt | awk '{print $1"\t"$2}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad.txt

    cat ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad.txt ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_strict.txt > ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad_strict.txt
    
    awk 'FNR==NR{seen[$2];next}($1 in seen){print}' ${PLINKDir}/${RunName}_ukb23157_c${chr}_v1_variant_filter.bim ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad_strict.txt > ${OutDir}/c${chr}_tmp && mv ${OutDir}/c${chr}_tmp ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad_strict.txt

    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_moderate.txt
    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_strict.txt
    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_cadd_pass_broad.txt
    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier.txt
    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_strict.txt
    rm ${OutDir}/${RunName}_c${chr}_b${block}_variants_impact_modifier_broad.txt

done

# Make a per-chromosome variant-gene file with annotation category

# Moderate - strict - canonical
cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_moderate_strict_canonical.txt | awk '{print $1"\t"$2"\tmoderate_cadd20_canonical"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt
cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_moderate_strict_alt_transcripts.txt | awk '{print $1"\t"$2"\tmoderate_cadd20_alt"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt
cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_moderate_broad_canonical.txt | awk '{print $1"\t"$2"\tmoderate_cadd1-20_canonical"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt
cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_moderate_broad_alt_transcripts.txt | awk '{print $1"\t"$2"\tmoderate_cadd1-20_alt"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_alt_transcripts.txt
cat ${OutDir}/${RunName}_c${chr}_b*_variants_impact_modifier_broad_strict.txt | awk '{print $1"\t"$2"\tmodifier_cadd1_all"}' | sort | uniq > ${OutDir}/${RunName}_c${chr}_variants_impact_modifier_broad_strict.txt

rm ${OutDir}/${RunName}_c${chr}_b*_variants_impact_*.txt

# Ensure all annotation categories are mutually exclusive (only one variant-gene combination). Prioritize the canonical and high-impact sets.

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt | awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' - ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt | awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' - ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt  > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt 

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt | awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' - ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt | awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' - ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_alt_transcripts.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_alt_transcripts.txt

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_alt_transcripts.txt | awk 'NR==FNR{a[$1"\t"$2];next} !($1"\t"$2 in a)' - ${OutDir}/${RunName}_c${chr}_variants_impact_modifier_broad_strict.txt > ${OutDir}/${RunName}_c${chr}_tmp && mv ${OutDir}/${RunName}_c${chr}_tmp ${OutDir}/${RunName}_c${chr}_variants_impact_modifier_broad_strict.txt


##################################################################################
## COMBINE SELECTED VARIANTS INTO FINAL ANNOTATION FILES FOR GENE-BASED TESTING ##
##################################################################################

cat ${OutDir}/${RunName}_c${chr}_variants_impact_high_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_high_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_strict_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_canonical.txt ${OutDir}/${RunName}_c${chr}_variants_impact_moderate_broad_alt_transcripts.txt ${OutDir}/${RunName}_c${chr}_variants_impact_modifier_broad_strict.txt > ${OutDir}/${RunName}_c${chr}_annotation_file.txt

awk '{print $1"\t"$2}' ${OutDir}/${RunName}_c${chr}_annotation_file.txt | sort | uniq > ${OutDir}/${RunName}_c${chr}_variant_gene_file.txt

for gene in $( awk '{print $2}' ${OutDir}/${RunName}_c${chr}_variant_gene_file.txt | sort | uniq ); do

    awk -v gene=${gene} '$2 == gene {print $1}' ${OutDir}/${RunName}_c${chr}_variant_gene_file.txt > ${OutDir}/${RunName}_c${chr}_tmp_variants.txt

    bp=$( sort ${OutDir}/${RunName}_c${chr}_tmp_variants.txt | head -1 | awk -F "_" '{print $2}' )

    echo -e ${gene}"\t"${chr}"\t"${bp} > ${OutDir}/${RunName}_c${chr}_set_list_tmp1.txt

    awk '
    { 
    for (i=1; i<=NF; i++)  {
    a[NR,i] = $i
    }
    }
    NF>p { p = NF }
    END {    
    for(j=1; j<=p; j++) {
    str=a[1,j]
    for(i=2; i<=NR; i++){
    str=str","a[i,j];
    }
    print str
    }
    }' ${OutDir}/${RunName}_c${chr}_tmp_variants.txt > ${OutDir}/${RunName}_c${chr}_set_list_tmp2.txt

    paste ${OutDir}/${RunName}_c${chr}_set_list_tmp1.txt ${OutDir}/${RunName}_c${chr}_set_list_tmp2.txt >> ${OutDir}/${RunName}_c${chr}_set_list.txt

    rm ${OutDir}/${RunName}_c${chr}_tmp_variants.txt
    rm ${OutDir}/${RunName}_c${chr}_set_list_tmp1.txt
    rm ${OutDir}/${RunName}_c${chr}_set_list_tmp2.txt

done

rm ${OutDir}/${RunName}_c${chr}_variant_gene_file.txt

