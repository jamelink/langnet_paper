cd /home/jitame/bin/software/mostest-master/

base_dir= "/data/clusterfs/lag/users/jitame/SENT_CORE";
apply_int=1;


mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_mean_fc")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_N29681_mean_fc_resid_mostest.txt") 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_mean_fc", "edges_chr"+i)
    mostest
    clear snps
end

clear nsubj
i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_mean_fc_resid_mostest.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_mean_fc", "edges_chr"+i)
mostest
clear snps

clear nsubj
i="XY"
pheno=fullfile(base_dir, "pheno", "sent_edges_N29681_mean_fc_resid_mostest_XY.txt");   
bfile=fullfile(base_dir, "geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir, "geno", "mostest", "out", "sent_edges_mean_fc", "edges_chr"+i)
mostest
clear snps



mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_asym_mean_hd")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_asym_N29681_mean_hd_resid_mostest.txt") 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_asym_mean_hd", "edges_asym_chr"+i)
    mostest
    clear snps
end

clear nsubj
i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_asym_N29681_mean_hd_resid_mostest.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_asym_mean_hd", "edges_asym_chr"+i)
mostest
clear snps

clear nsubj
i="XY"
pheno=fullfile(base_dir, "pheno", "sent_edges_asym_N29681_mean_hd_resid_mostest_XY.txt");   
bfile=fullfile(base_dir, "geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir, "geno", "mostest", "out", "sent_edges_asym_mean_hd", "edges_asym_chr"+i)
mostest
clear snps