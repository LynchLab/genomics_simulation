
#python rank_regions.py --bfile_train <bfile> --pheno_train <phenotype file> --meanLen <mean region length> --out <output file> --covar_train <covariates file>
python ~/src/MKLMM/rank_regions.py --bfile_train ../analysis_files/plink --pheno_train ../analysis_files/pheno --meanLen 1000 --out ../analysis_files/mklmm_regions.txt #--covar_train ../analysis_files/mklmm_cov.txt
python ~/src/MKLMM/mklmm_wrapper.py --bfile_train ../analysis_files/plink --pheno ../analysis_files/pheno --regions ../analysis_files/mklmm_regions.txt --numRegions 2 --kernel adapt --train_out ../analysis_files/predictions.txt
