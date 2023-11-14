chr=$1

out_dir="/data/clusterfs/lag/users/jitame/SENT_CORE/geno/regenie/exome/st2_out"

#mkdir -p ${out_dir}/c${chr}

cd ${out_dir}

dx cd /association_analyses/output/

dx download -r "1"

for gene in MANEAL TRIP11 DUSP29 DDX25 SLC25A48; do
dx download -r ${gene}
done
