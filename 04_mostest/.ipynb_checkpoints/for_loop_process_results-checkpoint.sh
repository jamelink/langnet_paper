pheno_path="/data/clusterfs/lag/users/jitame/SENT_CORE/pheno"
geno_path="/data/clusterfs/lag/users/jitame/SENT_CORE/geno/mostest/in"
mostest_path="/data/clusterfs/lag/users/jitame/SENT_CORE/geno/mostest/out"

python_path="/home/jitame/bin/envs/std_work_env/bin/python"


for i in {1..22}; do
    cd /home/jitame/bin/software/mostest-master/
    bim=$geno_path/mostest_geno_in_c${i}.bim
    mostest_res=$mostest_path/sent_edges_asym/edges_asym_chr${i}
    out=$mostest_res\_basic
    pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
    python process_results.py $bim $mostest_res $out
    #cd /home/jitame/bin/code/AICHA/05_mostest
    #out=$mostest_res/ext/
    #mkdir -p $out
    #$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out
done

for i in {1..22}; do
    cd /home/jitame/bin/software/mostest-master/
    bim=$geno_path/mostest_geno_in_c${i}.bim
    mostest_res=$mostest_path/sent_edges/edges_chr${i}
    out=$mostest_res\_basic
    pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
    python process_results.py $bim $mostest_res $out
    #cd /home/jitame/bin/code/AICHA/05_mostest
    #out=$mostest_res/ext/
    #mkdir -p $out
    #$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out
done

i=X
cd /home/jitame/bin/software/mostest-master/
bim=$geno_path/mostest_geno_in_c${i}.bim
mostest_res=$mostest_path/sent_edges_asym/edges_asym_chr${i}
out=$mostest_res\_basic
pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
python process_results.py $bim $mostest_res $out
#cd /home/jitame/bin/code/AICHA/05_mostest
#out=$mostest_res/ext/
#$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out

cd /home/jitame/bin/software/mostest-master/
bim=$geno_path/mostest_geno_in_c${i}.bim
mostest_res=$mostest_path/sent_edges/edges_chr${i}
out=$mostest_res\_basic
pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
python process_results.py $bim $mostest_res $out
#cd /home/jitame/bin/code/AICHA/05_mostest
#out=$mostest_res/ext/
#$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out

i=XY
cd /home/jitame/bin/software/mostest-master/
bim=$geno_path/mostest_geno_in_c${i}.bim
mostest_res=$mostest_path/sent_edges_asym/edges_asym_chr${i}
out=$mostest_res\_basic
pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
python process_results.py $bim $mostest_res $out
#cd /home/jitame/bin/code/AICHA/05_mostest
#out=$mostest_res/ext/
#$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out

cd /home/jitame/bin/software/mostest-master/
bim=$geno_path/mostest_geno_in_c${i}.bim
mostest_res=$mostest_path/sent_edges/edges_chr${i}
out=$mostest_res\_basic
pheno=$pheno_path/sent_edges_N29681_resid_mostest.txt
python process_results.py $bim $mostest_res $out
#cd /home/jitame/bin/code/AICHA/05_mostest
#out=$mostest_res/ext/
#$python_path 04_process_results_ext_JS.py $bim $mostest_res $pheno $out


cat $mostest_path/sent_edges_asym/edges_asym_chr*_basic.most_orig.sumstats > $mostest_path/sent_edges_asym/edges_asym.sumstats
cat $mostest_path/sent_edges/edges_chr*_basic.most_orig.sumstats > $mostest_path/sent_edges/edges.sumstats