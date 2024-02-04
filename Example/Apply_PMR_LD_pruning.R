#!/usr/bin/env Rscript


Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

#library(gwplotting)
library(PMR)
library(data.table)
library(readr)
library(dplyr)
library(SumTool)




################################
# 1. read in pQTL data

pQTL = fread(pQTL_per_chr)
pQTL = pQTL %>% filter(TargetID == gene_name)

pQTL$ID = paste(pQTL$CHROM,pQTL$POS,pQTL$A1,pQTL$A2,sep = ":")
pQTL$IDflip = paste(pQTL$CHROM,pQTL$POS,pQTL$A2,pQTL$A1,sep = ":")

################################
# 2. read in GWAS data
gwas <- fread(gwas_file,header = T)

names(gwas)[1] <- "CHROM"
gwas$ID <- paste(gwas$CHROM,gwas$POS,gwas$REF,gwas$ALT,sep=":")
gwas$IDflip <- paste(gwas$CHROM,gwas$POS,gwas$ALT,gwas$REF,sep=":")


###########################################
# 3. intersect pQTL-geno-GWAS SNPs

# read in genotype data 
geno <- read_binary(bfile=geno_file, additive=TRUE, threads=1)

ref.geno <- geno$geno
ref.map <- geno$map
# create geno SNP ID
ref.map$SNP = paste(ref.map$Chr,ref.map$Pos,ref.map$Ref,ref.map$Alt,sep=":")

# intersect snps of pQTL-geno
pQTL_snp_noflip <- intersect(ref.map$SNP,pQTL$ID) 
pQTL_snp_flip <- intersect(ref.map$SNP,pQTL$IDflip) 

pQTL_inter_snp <- c(pQTL_snp_noflip,pQTL_snp_flip)
pQTL_snp_flip_index <- which(pQTL$IDflip  %in% pQTL_snp_flip) 

# intersect snps of GWAS-geno
gwas_snp_noflip <-  intersect(ref.map$SNP,gwas$ID) 
gwas_snp_flip <- intersect(ref.map$SNP,gwas$IDflip) 

gwas_inter_snp <- c(gwas_snp_noflip,gwas_snp_flip)
gwas_snp_flip_index <- which(gwas$IDflip  %in% gwas_snp_flip)

# intersect of pQTL-GWAS-geno
inter_snp <- intersect(pQTL_inter_snp,gwas_inter_snp) 



###########################################
# 4. filter SNPs in genotype, LD pruning calculate LD matrix

# filter snps in genotype
ref.map.filter <- ref.map %>% filter(SNP %in% inter_snp)

geno_snp_index1 <- which(ref.map$SNP %in% inter_snp)
ref.geno.filter <- deepcopy(ref.geno, cols=geno_snp_index1)


# do the LD pruning and get the index
snp_prune <- LDprune(geno = ref.geno.filter, map = ref.map.filter,
                     r2.cutoff = 0.8, w = 100000, threads = 1)

geno_snp_index2 <- which(ref.map.filter$SNP %in% snp_prune)
ref.geno.prune <- deepcopy(ref.geno.filter, cols=geno_snp_index2)



# caculate LD matrix1
cat("calculate LD matrix\n")

LD <- LDcal(geno= ref.geno.prune, threads=1)


###########################################
# 5. filter LD pruning SNPs in pQTL and gwas

# replace the ID of flipped pQTLs
pQTL$ID[pQTL_snp_flip_index] = pQTL$IDflip[pQTL_snp_flip_index] 

# flip the sigh of T-stat and BETA for flipped snps
pQTL$Zscore[pQTL_snp_flip_index] = -1*pQTL$Zscore[pQTL_snp_flip_index]


# calculate Zscore from T-stat

# filter LD pruning pQTL snps
pQTL_filter = pQTL %>% filter(ID %in% snp_prune)


###
# replace the ID of flipped gwas SNPs
gwas$ID[gwas_snp_flip_index] = gwas$IDflip[gwas_snp_flip_index] 

# flip the sigh of T-stat and BETA for flipped snps
gwas$Zscore[gwas_snp_flip_index] = -1*gwas$Zscore[gwas_snp_flip_index]

# filter LD pruning gwas snps
gwas_filter = gwas %>% filter(ID %in% snp_prune)

# reaplce inf value in gwas data
gwas_inf_index <- which(is.infinite(gwas_filter$Zscore))
#gwas_Zscore_vec[gwas_inf_index] <- max(gwas_Zscore_vec[-gwas_inf_index])
gwas_filter$Zscore[which(gwas_filter$Zscore==Inf)] <- max(gwas_filter$Zscore[-gwas_inf_index])
gwas_filter$Zscore[which(gwas_filter$Zscore==-Inf)] <- min(gwas_filter$Zscore[-gwas_inf_index])


########################################
# 6. specify sample size for pQTL and gwas data
# load the sample size n1 from eQTL data and n2 from GWAS data 
n1= mean(pQTL_filter$N) 
n2= floor(mean(gwas_filter$N))   

###################################
# 7. conduct PMR
pQTL_Zscore_vec = pQTL_filter$Zscore
gwas_Zscore_vec = gwas_filter$Zscore

cat("conducting summary level PMR\n")
result = PMR_summary_Egger(pQTL_Zscore_vec, gwas_Zscore_vec, LD, LD, 
                           samplen1=n1, samplen2=n2, 
                           lambda=0, max_iterin =1000,epsin=1e-5, 
                           Heritability_geneexpression_threshold=0)

# save result file
result$chr = chr
result$gene = gene_name
result$nsnp_pQTL= length(pQTL_Zscore_vec)
result$nsnp_gwas = length(gwas_Zscore_vec)
result$gwas_N = n2


#out_dir="/projects/YangLabData/thu/pQTL_sum/PMR_pQTL/output/brain"
cat(paste0("saving result into"),out_dir,"\n")
saveRDS(result, paste0(out_dir,"/",gene_name,"_PMR_result_LD_0.8.rds"))
