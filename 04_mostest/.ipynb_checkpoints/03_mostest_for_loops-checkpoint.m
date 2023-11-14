cd /home/jitame/bin/software/mostest-master/

base_dir="/data/clusterfs/lag/users/jitame/SENT_CORE";
apply_int=1;

mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_L")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_L_mostest.txt"); 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out","sent_edges_L", "edges_L_chr"+i)
    mostest
    clear snps
end

mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_R")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_R_mostest.txt") 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_R", "edges_R_chr"+i)
    mostest
    clear snps
end

mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_intra")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_intra_mostest.txt"); 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out","sent_edges_intra", "edges_intra_chr"+i)
    mostest
    clear snps
end

mkdir fullfile(base_dir, "geno", "mostest", "out", "sent_edges_inter")
for i=1:22
    pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_inter_mostest.txt") 
    bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
    out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_inter", "edges_inter_chr"+i)
    mostest
    clear snps
end

i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_intra_mostest.txt"); 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out","sent_edges_intra", "edges_intra_chr"+i)
mostest
clear snps

i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_inter_mostest.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_inter", "edges_inter_chr"+i)
mostest
clear snps

i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_R_mostest.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_R", "edges_R_chr"+i)
mostest
clear snps

i="X"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_L_mostest.txt"); 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out","sent_edges_L", "edges_L_chr"+i)
mostest
clear snps

clear nsubj

i="XY"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_intra_mostest_XY.txt"); 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out","sent_edges_intra", "edges_intra_chr"+i)
mostest
clear snps

i="XY"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_inter_mostest_XY.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_inter", "edges_inter_chr"+i)
mostest
clear snps

i="XY"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_R_mostest_XY.txt") 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out", "sent_edges_R", "edges_R_chr"+i)
mostest
clear snps

i="XY"
pheno=fullfile(base_dir,"pheno","sent_edges_N29681_resid_heritable_edges_L_mostest_XY.txt"); 
bfile=fullfile(base_dir,"geno", "mostest", "in", "mostest_geno_in_c"+i)
out=fullfile(base_dir,"geno","mostest", "out","sent_edges_L", "edges_L_chr"+i)
mostest
clear snps

