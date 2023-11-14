base_dir=/data/clusterfs/lag/users/jitame/SENT_CORE

out_path=/data/clusterfs/lag/users/jitame/SENT_CORE/geno/polygenic-scores/prs_estim/


for input_ma in $(ls /data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/gwas_summary/all_sums_pgs/*.ma); do
pheno_name=$(basename ${input_ma})
pheno_name=${pheno_name::-3}
echo $pheno_name

#mkdir $out_path/$pheno_name
#done

for i in {1..22};
do
qsub -q multi15.q -N "prs_cs_${pheno_name}_${i}" /home/jitame/bin/code/AICHA/07_pgs/prs_1.sh -f /data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/gwas_summary/all_sums_pgs/$pheno_name.ma -n ${cat /data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/gwas_summary/all_sums_pgs/${pheno_name}_n.txt} -o /data/clusterfs/lag/users/jitame/SENT_CORE/geno/polygenic-scores/prs_estim/${pheno_name}/${pheno_name} -c $i 
done

done



## ================================================================================== ##

#PRS 2

base_dir=/data/clusterfs/lag/users/jitame/SENT_CORE

for input in asd dyslexia hand read scz; do
mkdir $base_dir/geno/polygenic-scores/prs_out/${input}

for i in {1..22}; 
do
qsub -q single15.q -N "prs_plink_${input}_${i}" /home/jitame/bin/code/AICHA/07_pgs/prs_2.sh -f $base_dir/geno/polygenic-scores/prs_estim/${input}/${input}_pst_eff_a1_b0.5_phi1e-02_chr${i}.txt  -n $base_dir/geno/polygenic-scores/prs_out/${input}/${input}_prs_chr${i} -c $i 
done

done
