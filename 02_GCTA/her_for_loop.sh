#!/bin/sh

base_dir=/data/clusterfs/lag/users/jitame/SENT_CORE
wrapper=/home/jitame/bin/code/AICHA/02_GCTA/01_run_heritability.sh

mkdir -p $base_dir/geno/h2/edges
for i in {1..630}; do
  out_hsq="$base_dir/geno/h2/edges/sent_edges_${i}.hsq"
  if [ -f ${out_hsq} ]; then
  echo "${out_hsq} exists"
  else
  echo "${out_hsq} does not exist"
qsub -q multi15.q -p -50 -N h2_s_e_${i} $wrapper -f $base_dir/pheno/sent_edges_N29681_resid_norm.txt -n $i -c $base_dir/geno/h2/edges/sent_edges_${i}
  fi
done

mkdir -p $base_dir/geno/h2/edges_asym
for i in {1..153}; do
  out_hsq="$base_dir/geno/h2/edges_asym/sent_edges_asym_${i}.hsq"
  if [ -f ${out_hsq} ]; then
  echo "${out_hsq} exists"
  else
  echo "${out_hsq} does not exist"
  qsub -q multi15.q -p -50 -N h2_s_e_asym_${i} $wrapper -f $base_dir/pheno/sent_edges_asym_N29681_resid_norm.txt -n $i -c $base_dir/geno/h2/edges_asym/sent_edges_asym_${i}
 fi
done
