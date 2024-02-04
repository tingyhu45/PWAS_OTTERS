# PWAS_OTTERS
Proteome-wide association study (PWAS) leveraging pQTL summary data

### Getting started

* Download [OTTERS](https://github.com/daiqile96/OTTERS)  using:
  ```
  git clone https://github.com/daiqile96/OTTERS.git
  ```

* Download the required software and packages
  * [BGZIP](http://www.htslib.org/doc/bgzip.html)
  * [TABIX](http://www.htslib.org/doc/tabix.html)
  * [PLINK 1.9](https://www.cog-genomics.org/plink/)
  * Python modules/libraries: pandas 1.4.4, scipy 1.7.3, numpy 1.21.5, pysam 0.19.1


* To apply SDPR and lassosum as imputation models
  * [SDPR](https://github.com/eldronzhou/SDPR)
  * R packages
    * [lassosum](https://github.com/tshmak/lassosum) to perform lassosum
    * [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html) to perform pseudovalidation in lassosum

Detailed instructions for installment can be found at [OTTERS](https://github.com/daiqile96/OTTERS/blob/main/README.md)

* Download R pacakge [PMR-Egger](https://github.com/yuanzhongshang/PMR) for validation of PWAS risk genes
  

### Example

* Apply OTTERS
  ```
  # activate the environment
  conda activate otters
  
  # set number of threads to be used
  N_THREADS=1
  
  # set up my OTTERS directory and SDPR directory
  OTTERS_DIR=/home/qdai8/projects/bin/OTTERS
  SDPR_DIR=/home/qdai8/projects/bin/SDPR
  
  # make sure the dynamic libraries of SDPR are not changed (For SDPR)
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
  
  # prevent automatically using  all available cores on a compute node (For SDPR and PRS-CS)
  export MKL_NUM_THREADS=$N_THREADS
  export NUMEXPR_NUM_THREADS=$N_THREADS
  export OMP_NUM_THREADS=$N_THREADS
  
  # Load R to perform lassosum
  module load R
  
  # Start to run OTTERS
  cd ${OTTERS_DIR}/Example
  
  # Input for OTTERS STAGE I 
  # Annotation File 
  exp_anno=exp_anno.txt
  # Genotype data from LD reference panel
  geno_dir=Exp_geno
  # eQTL summary statistics 
  sst_file=Exp_eQTLSumStats.txt
  # Input for OTTERS STAGE II
  # GWAS summary statistics 
  gwas_sst_file=Exp_GWASSumStats.txt
  
  # Set chromosome number (The example is for chromosome 4)
  chr=4
  # Set LD-clumping threshold in STAGE I
  clump_r2=0.99
  # Set output directory for STAGE I
  out_dir=Results
  
  # STAGE I
  # train eQTL weights using P+T, lassosum, SDPR and PRS-CS. 
  # It may take several minutes to complete.
  python3 ${OTTERS_DIR}/training.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --SDPR_dir=${SDPR_DIR} \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --sst_file=${sst_file} \
  --out_dir=${out_dir} \
  --chrom=${chr} \
  --r2=${clump_r2} \
  --models=PT,lassosum,SDPR,PRScs \
  --thread=$N_THREADS
  
  # Set output directory for STAGE II
  twas_dir=TWAS
  
  # STAGE II
  # gene-based association test using eQTL-weight trained from P+T, lassosum, SDPR and PRS-CS.
  python3 ${OTTERS_DIR}/testing.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --weight_dir=${OTTERS_DIR}/Example/Results \
  --models=P0.001,P0.05,lassosum,SDPR,PRScs \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --out_dir=${twas_dir} \
  --gwas_file=${gwas_sst_file} \
  --chrom=${chr} \
  --thread=$N_THREADS
  
  # get imputed genetically regulated gene expression
  impute_dir=GReX
  # samples to perform imputation
  samples=Exp_samples.txt
  # imputation
  python3 ${OTTERS_DIR}/imputing.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --weight_dir=${OTTERS_DIR}/Example/Results \
  --models=P0.001,P0.05,lassosum,SDPR,PRScs \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --out_dir=${impute_dir} \
  --chrom=${chr} \
  --samples=${samples} \
  --thread=$N_THREADS
  ```

  
* Combine p-values using ACAT-O test
Launch R and load the package
```
  library(ACAT)
  
  ACAT_withNA = function(p_vec){
    p_vec_noNA = p_vec[is.na(p_vec) == F]
    ACAT(p_vec_noNA)
    
  }
  ```
  
* Apply PMR-Egger for validating the causality of PWAS risk genes
  Launch R and load the package
  ```
  library(PMR)
  
  PMR_summary_Egger(pQTL_Zscore_vec, gwas_Zscore_vec, LD, LD, 
                           samplen1=n1, samplen2=n2, 
                           lambda=0, max_iterin =1000,epsin=1e-5, 
                           Heritability_geneexpression_threshold=0)
  ```
  * The code for applying PMR-Egger with LD clumping can be found at 
