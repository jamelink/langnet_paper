#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -m eas
#$ -N subset_snpstats
#$ -q single15.q


for c in {1..22}; do
qsub -N "subset_40k_c${c}" /home/jitame/bin/code/AICHA/genetics/regenie/gwas/imaging40k_subset_and_snpstats.sh \
-c ${c} \
-s /home/jitame/bin/code/AICHA/genetics/regenie/gwas/imaging40k_subset_and_snpstats_config.txt 
done

for c in {1..18}; do
qsub -N "variant_qc_40k_c${c}" /home/jitame/bin/code/AICHA/genetics/regenie/gwas/variant_qc_wrapper.sh \
-c ${c}
done

for c in {1..22}; do
qsub -N "mkbed_c${c}_mostest" /home/jitame/bin/code/AICHA/genetics/mostest/convert_bgen_mostest.sh \
-c ${c}
done