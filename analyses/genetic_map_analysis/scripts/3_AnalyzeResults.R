#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_AnalyzeResults.R
# ------------------
# Analyze haplotype block results.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("cowplot")
library("mltools")
library("qqman")
library("ggpubr")
library("ggvenn")
library("RColorBrewer")
library("plotROC")
library("ggcorrplot")
library("corrplot")
library("gridExtra")
library("igraph")
library("cluster")
library("factoextra")
library("pvclust")

### Mean wrapper function ###
my_mean <- function( x )
{
  return(mean(x, na.rm=T))
}

### Sum wrapper function ###
my_sum <- function( x )
{
  return(sum(x, na.rm=T))
}

### Minimum wrapper function ###
my_min <- function( x )
{
  return(min(x, na.rm=T))
}

### Maximum wrapper function ###
my_max <- function( x )
{
  return(max(x, na.rm=T))
}

### Load eQTLs ###
load_eQTLs <- function( phenotype )
{
  #--------------------------#
  # 1) Load eQTLs data       #
  #--------------------------#
  EQTLs    = readRDS(paste0("./data/tribolium_eqtl/significant/HD_G1_Tcas3.30_imputed_",phenotype,"_significant.rds"))
  EQTLs$ID = EQTLs$rs
  #--------------------------#
  # 2) Load SNP annotation   #
  #--------------------------#
  ANNOTATION        = read.table("./data/tribolium_snp/snp_table_ALL_Tcas3.30_raw_SNP.csv", h=T, sep="\t")
  ANNOTATION$gene   = ANNOTATION$Feature_id
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "Annotation", "Putative_impact", "gene"))
  #--------------------------#
  # 3) Select and merge data #
  #--------------------------#
  eQTLs = merge(EQTLs, ANNOTATION_SELECT, by="ID")
  #--------------------------#
  # 4) Return data           #
  #--------------------------#
  rm(ANNOTATION)
  return(eQTLs)
}

### Apply genomic control adjustment on p-values ###
calculate_gc <- function( p )
{
  n      = length(p)
  x2obs  = qchisq(p, 1, lower.tail=FALSE)
  x2exp  = qchisq(1:n/n, 1, lower.tail=FALSE)
  lambdA = median(x2obs)/median(x2exp)
  x2new  = x2obs/lambdA
  gc     = pchisq(x2new, df=1, lower.tail=FALSE)
  return(gc)
}

### Calculate the FDR ###
calculate_fdr <- function( p )
{
  fdr = p.adjust(p, method="fdr")
  return(fdr)
}

### Load and merge permutation results ###
load_permutation_results <- function( HB_dataset, classification_score )
{
  #---------------------------------------#
  # 1) Load and merge permutation results #
  #---------------------------------------#
  VARIABLES = c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift",
                "evol_response_pool", "eQTL_carrier", "eQTL_phenotype", "is_ASE", "significant_sg_pool",
                "hub_gene", "highest_sg_pool", "highest_de_pool", "abs_logFC", "log_net_selection")
  MERGED    = data.frame()
  first     = T
  for (variable in VARIABLES)
  {
    d = read.table(paste0("./analyses/genetic_map_analysis/data/permutation_tests/",variable,"_merged.txt"), sep=";", h=T)
    if (first)
    {
      MERGED        = data.frame(d[,c("index", "haplotype_block", "pvalue")])
      MERGED$gc     = calculate_gc(MERGED$pvalue)
      MERGED$fdr    = calculate_fdr(MERGED$pvalue)
      MERGED$fdr_gc = calculate_fdr(MERGED$gc)
      if (classification_score == "pvalue")
      {
        MERGED$signif = as.numeric(MERGED$pvalue < 0.05)
      } else if (classification_score == "bonferroni") {
        N             = length(MERGED$pvalue)
        MERGED$signif = as.numeric(MERGED$pvalue < 0.05/N)
      } else if (classification_score == "gc") {
        MERGED$signif = as.numeric(MERGED$gc < 0.05) 
      } else if (classification_score == "fdr") {
        MERGED$signif = as.numeric(MERGED$fdr < 0.05)
      } else if (classification_score == "fdr_gc") {
        MERGED$signif = as.numeric(MERGED$fdr_gc < 0.05)
      }
      names(MERGED) = c("index", "haplotype_block", paste0(variable,"_pvalue"), paste0(variable,"_gc"), paste0(variable,"_fdr"), paste0(variable,"_fdr_gc"), variable)
      first         = F
    } else {
      d        = d[,c("index", "pvalue")]
      d$gc     = calculate_gc(d$pvalue)
      d$fdr    = calculate_fdr(d$pvalue)
      d$fdr_gc = calculate_fdr(d$gc)
      if (classification_score == "pvalue")
      {
        d$signif = as.numeric(d$pvalue < 0.05)
      } else if (classification_score == "bonferroni") {
        N        = length(d$pvalue)
        d$signif = as.numeric(d$pvalue < 0.05/N)
      } else if (classification_score == "gc") {
        d$signif = as.numeric(d$gc < 0.05)
      } else if (classification_score == "fdr") {
        d$signif = as.numeric(d$fdr < 0.05)
      } else if (classification_score == "fdr_gc") {
        d$signif = as.numeric(d$fdr_gc < 0.05)
      }
      names(d) = c("index", paste0(variable,"_pvalue"), paste0(variable,"_gc"), paste0(variable,"_fdr"), paste0(variable,"_fdr_gc"), variable)
      MERGED   = merge(MERGED, d, by="index")
    }
  }
  #---------------------------------------#
  # 2) Compute parallelism variables      #
  #---------------------------------------#
  MERGED$Parallelism                = MERGED$L1_shift+MERGED$L2_shift+MERGED$L3_shift+MERGED$L5_shift+MERGED$L6_shift+MERGED$Mx1_shift+MERGED$Mx2_shift
  MERGED$NoShift                    = as.numeric(MERGED$Parallelism == 0)
  MERGED$isShifting                 = as.numeric(MERGED$Parallelism >= 1)
  MERGED$isSingleShifting           = as.numeric(MERGED$Parallelism == 1)
  MERGED$isPartiallyParallel        = as.numeric(MERGED$Parallelism > 1 & MERGED$Parallelism < 5)
  MERGED$isPartiallyorFullyParallel = as.numeric(MERGED$Parallelism > 1)
  MERGED$isFullyParallel            = as.numeric(MERGED$Parallelism > 4)
  MERGED$P0                         = as.numeric(MERGED$Parallelism == 0)
  MERGED$P1                         = as.numeric(MERGED$Parallelism == 1)
  MERGED$P2                         = as.numeric(MERGED$Parallelism == 2)
  MERGED$P3                         = as.numeric(MERGED$Parallelism == 3)
  MERGED$P4                         = as.numeric(MERGED$Parallelism == 4)
  MERGED$P5                         = as.numeric(MERGED$Parallelism == 5)
  MERGED$P6                         = as.numeric(MERGED$Parallelism == 6)
  MERGED$P7                         = as.numeric(MERGED$Parallelism == 7)
  #---------------------------------------#
  # 3) Collect additional information     #
  #---------------------------------------#
  ### 3.1) Chromosome location ###
  hb_chr            = HB_dataset[!duplicated(HB_dataset$haplotype_block),c("haplotype_block", "CHROM")]
  row.names(hb_chr) = hb_chr$haplotype_block
  MERGED$chr        = hb_chr[as.character(MERGED$haplotype_block), "CHROM"]
  MERGED            = MERGED[order(MERGED$chr, MERGED$haplotype_block),]
  ### 3.2) Custom index per chromosome ###
  current_chr = MERGED$chr[1]
  counter     = 1
  chr_index   = c()
  for(i in seq(1, dim(MERGED)[1]))
  {
    new_chr = MERGED$chr[i]
    if (current_chr == new_chr)
    {
      chr_index = c(chr_index, counter)
      counter   = counter+1
    } else {
      current_chr = new_chr
      counter     = 1
      chr_index = c(chr_index, counter)
    }
  }
  MERGED$chr_index = chr_index
  MERGED$chr_int   = as.numeric(as.factor(MERGED$chr))
  ### 3.3) Mean |log(FC)| ###
  mean_logFC        = tapply(HB_dataset$abs_logFC, HB_dataset$haplotype_block, my_mean)
  MERGED$mean_logFC = mean_logFC[as.character(MERGED$haplotype_block)]
  ### 3.4) Max |log(FC)| ###
  max_logFC        = tapply(HB_dataset$abs_logFC, HB_dataset$haplotype_block, my_max)
  MERGED$max_logFC = max_logFC[as.character(MERGED$haplotype_block)]
  ### 3.5) Mean log(|Sg|) ###
  mean_sg        = tapply(HB_dataset$log_net_selection, HB_dataset$haplotype_block, my_mean)
  MERGED$mean_sg = mean_sg[as.character(MERGED$haplotype_block)]
  ### 3.6) Max log(|Sg|) ###
  max_sg        = tapply(HB_dataset$log_net_selection, HB_dataset$haplotype_block, my_max)
  MERGED$max_sg = max_sg[as.character(MERGED$haplotype_block)]
  #---------------------------------------#
  # 4) Clean the dataset                  #
  #---------------------------------------#
  MERGED                    = MERGED[,c("haplotype_block", "L1_shift", "L2_shift", "L3_shift", "L5_shift",
                                        "L6_shift", "Mx1_shift", "Mx2_shift",
                                        "evol_response_pool", "highest_de_pool", "eQTL_carrier", "eQTL_phenotype",
                                        "abs_logFC", "log_net_selection", "is_ASE",
                                        "significant_sg_pool", "highest_sg_pool", "hub_gene",
                                        "Parallelism", "NoShift", "isShifting", "isSingleShifting", "isPartiallyParallel",
                                        "isPartiallyorFullyParallel", "isFullyParallel",
                                        "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7")]
  rownames(MERGED)          = MERGED$haplotype_block
  to_correct                = which(names(MERGED)=="evol_response_pool")
  names(MERGED)[to_correct] = "significant_de_pool"
  to_correct                = which(names(MERGED)=="abs_logFC")
  names(MERGED)[to_correct] = "p_abs_logFC"
  to_correct                = which(names(MERGED)=="log_net_selection")
  names(MERGED)[to_correct] = "p_log_net_selection"
  #---------------------------------------#
  # 5) Return the final merged data       #
  #---------------------------------------#
  return(MERGED)
}

### Merge haplotype blocks and permutation tests datasets ###
merge_HB_permutations_datasets <- function( HB_dataset, PERM_dataset )
{
  to_keep        = c("haplotype_block", "hb_size", "ID", "POS", "CHROM", "gene", "logFC", "abs_logFC", "net_selection", "direct_selection", "diff_selection", "log_net_selection", "log_direct_selection")
  HB_dataset     = HB_dataset[,to_keep]
  MERGED_dataset = merge(HB_dataset, PERM_dataset, by="haplotype_block")
  return(MERGED_dataset)
}

### Merge haplotype block and eQTLs datasets ###
merge_HB_eQTLs_datasets <- function( HB_PERM_dataset, eQTLs )
{
  Pleiotropy          = c()
  Connectivity        = c()
  isPleiotropic       = c()
  isHighlyPleiotropic = c()
  isMultiRegulated    = c()
  isHighlyRegulated   = c()
  for (hb in unique(HB_PERM_dataset$haplotype_block))
  {
    dl       = filter(HB_PERM_dataset, haplotype_block==hb)
    N        = dim(dl)[1]
    carrier  = unique(dl$eQTL_carrier)
    pheno    = unique(dl$eQTL_phenotype)
    pleio    = 0
    connec   = 0
    ispleio  = 0
    isreg    = 0
    isHpleio = 0
    isHreg   = 0
    if (carrier == 1)
    {
      phenotypes = unique(filter(eQTLs, ID %in% dl$ID)$phenotype)
      pleio      = length(unique(filter(HB_PERM_dataset, gene%in%phenotypes)$haplotype_block))
      ispleio    = as.numeric(pleio > 1)
      isHpleio   = as.numeric(pleio > 4)
    }
    if (pheno == 1)
    {
      genes  = unique(filter(eQTLs, phenotype %in% dl$gene)$gene)
      connec = length(unique(filter(HB_PERM_dataset, gene%in%genes)$haplotype_block))
      isreg  = as.numeric(connec > 1)
      isHreg = as.numeric(connec > 4)
    }
    Pleiotropy          = c(Pleiotropy, rep(pleio, N))
    Connectivity        = c(Connectivity, rep(connec, N))
    isPleiotropic       = c(isPleiotropic, rep(ispleio, N))
    isHighlyPleiotropic = c(isHighlyPleiotropic, rep(isHpleio, N))
    isMultiRegulated    = c(isMultiRegulated, rep(isreg, N))
    isHighlyRegulated   = c(isHighlyRegulated, rep(isHreg, N))
  }
  HB_PERM_dataset$Pleiotropy          = Pleiotropy
  HB_PERM_dataset$Connectivity        = Connectivity
  HB_PERM_dataset$isPleiotropic       = isPleiotropic
  HB_PERM_dataset$isHighlyPleiotropic = isHighlyPleiotropic
  HB_PERM_dataset$isMultiRegulated    = isMultiRegulated
  HB_PERM_dataset$isHighlyRegulated   = isHighlyRegulated
  return(HB_PERM_dataset)
}

### Run matrix tests (for now, Exact Fisher tests) ###
run_matrix_tests <- function( dataset, variables_to_test, chromosome )
{
  D = data.frame(dataset)
  if (chromosome > 0)
  {
    D = filter(D, chr_int==chromosome)
  }
  N                = length(variables_to_test)
  Mpearson         = matrix(rep(1.0, N*N), ncol=N)
  Mpearson_signif  = matrix(rep(0.0, N*N), ncol=N)
  Mspearman        = matrix(rep(1.0, N*N), ncol=N)
  Mspearman_signif = matrix(rep(0.0, N*N), ncol=N)
  Mfisher          = matrix(rep(1.0, N*N), ncol=N)
  Mfisher_signif   = matrix(rep(0.0, N*N), ncol=N)
  for(i in seq(1,N))
  {
    for(j in seq(1,N))
    {
      var1                  = variables_to_test[i]
      var2                  = variables_to_test[j]
      cor_coef              = cor.test(D[,var1], D[,var2])$estimate
      pearson_pval          = cor.test(D[,var1], D[,var2], method="pearson")$p.value
      spearman_pval         = cor.test(D[,var1], D[,var2], method="spearman")$p.value
      Mpearson[i,j]         = pearson_pval
      Mpearson_signif[i,j]  = as.numeric(pearson_pval<0.05)
      Mspearman[i,j]        = spearman_pval
      Mspearman_signif[i,j] = as.numeric(spearman_pval<0.05)
      test_table            = table(D[,var1], D[,var2])
      if (!(dim(test_table)[1]==2 & dim(test_table)[2]==2))
      {
        Mfisher[i,j]        = 1.0
        Mfisher_signif[i,j] = 0.0
      } else {
        fisher_pval         = fisher.test(test_table)$p.value
        Mfisher[i,j]        = fisher_pval
        Mfisher_signif[i,j] = as.numeric(fisher_pval<0.05)
      }
      if (!is.na(cor_coef) & cor_coef < 0.0)
      {
        #Mpearson[i,j]         = -Mpearson[i,j]
        Mpearson_signif[i,j]  = -Mpearson_signif[i,j]
        #Mspearman[i,j]        = -Mspearman[i,j]
        Mspearman_signif[i,j] = -Mspearman_signif[i,j]
        #Mfisher[i,j]          = -Mfisher[i,j]
        Mfisher_signif[i,j]   = -Mfisher_signif[i,j]
      }
    }
    colnames(Mpearson)         = variables_to_test
    rownames(Mpearson)         = variables_to_test
    colnames(Mpearson_signif)  = variables_to_test
    rownames(Mpearson_signif)  = variables_to_test
    colnames(Mspearman)        = variables_to_test
    rownames(Mspearman)        = variables_to_test
    colnames(Mspearman_signif) = variables_to_test
    rownames(Mspearman_signif) = variables_to_test
    colnames(Mfisher)          = variables_to_test
    rownames(Mfisher)          = variables_to_test
    colnames(Mfisher_signif)   = variables_to_test
    rownames(Mfisher_signif)   = variables_to_test
  }
  return(list("pearson"=Mpearson, "pearson_signif"=Mpearson_signif,
              "spearman"=Mspearman, "spearman_signif"=Mspearman_signif,
              "fisher"=Mfisher, "fisher_signif"=Mfisher_signif))
}

### Aggregator function ###
aggregator_function <- function( x, aggregator )
{
  if (aggregator == "sum")
  {
    return(sum(x, na.rm=T))
  }
  if (aggregator == "mean")
  {
    return(mean(x, na.rm=T))
  }
  if (aggregator == "median")
  {
    return(median(x, na.rm=T))
  }
}

### Run one permutations test ###
run_permutations_test <- function( dataset, target_variable, source_variable, source_list, nb_reps, aggregator )
{
  x              = dataset[,target_variable]
  pos            = which(dataset[,source_variable]%in%source_list)
  original_score = aggregator_function(x[pos], aggregator)
  permutations   = t(sapply(1:nb_reps,function(i){
    score = aggregator_function(sample(x)[pos], aggregator)
    return(score)
  }))
  qval = sum(permutations<original_score)/nb_reps
  return(1-qval)
}



##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")


#----------------------------#
# 1) Load and build datasets #
#----------------------------#
#eQTLs           = load_eQTLs("EXPRESSION")
HB_dataset      = readRDS(file="./analyses/genetic_map_analysis/data/haplotype_block_dataset.rds")
PERM_dataset    = load_permutation_results(HB_dataset, "fdr")
HB_PERM_dataset = merge_HB_permutations_datasets(HB_dataset, PERM_dataset)
MERGED_dataset  = merge_HB_eQTLs_datasets(HB_PERM_dataset, eQTLs)
MERGED_dataset  = MERGED_dataset[!duplicated(MERGED_dataset$haplotype_block),]

MERGED_dataset = filter(MERGED_dataset, hb_size >= 10)

#----------------------------#
# 2) Compute matrix tests    #
#----------------------------#
TO_TEST = c("NoShift", "isShifting", "isPartiallyorFullyParallel", "isFullyParallel",
            "significant_de_pool", "highest_de_pool", "p_abs_logFC",
            "eQTL_carrier", "eQTL_phenotype", "isPleiotropic", "isMultiRegulated",
            "isHighlyPleiotropic", "isHighlyRegulated")
RES     = run_matrix_tests(MERGED_dataset, TO_TEST, 0)
ggcorrplot(RES$fisher_signif, type="upper")

#----------------------------#
# 3) Look at shifting blocks #
#----------------------------#
M           = MERGED_dataset[,c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift")]
res.dist    = get_dist(M, method="binary")
res.hc      = hclust(res.dist)
#plot(res.hc, cex = 0.5)
clusts      = cutree(res.hc, h=0.1)
BLOCK_DATA  = data.frame(MERGED_dataset)
pheno_var   = c("evol_response_pool", "eQTL_carrier", "eQTL_phenotype", "is_ASE", "hub_gene")
TEST_RES    = data.frame()
Cluster     = c()
Parallelism = c()
for(clust in sort(unique(clusts)))
{
  print(clust)
  indices   = names(clusts[clusts==clust])
  hbs       = filter(BLOCK_DATA, index%in%indices)$haplotype_block
  col_name              = paste0("clust_",clust)
  BLOCK_DATA[,col_name] = as.numeric(BLOCK_DATA$index%in%indices)
  TO_TEST     = c(col_name, pheno_var)
  RES         = run_matrix_tests(BLOCK_DATA, TO_TEST, 0)
  TEST_RES    = rbind(TEST_RES, RES$fisher_signif[col_name,-1])
  Cluster     = c(Cluster, clust)
  Parallelism = c(Parallelism, mean(filter(BLOCK_DATA, index%in%indices)$Parallelism))
  # nb_eQTL = sum(filter(BLOCK_DATA, index%in%indices)$eQTL_carrier)
  # nb_DE = sum(filter(BLOCK_DATA, index%in%indices)$evol_response_pool)
  # nb_HUB = sum(filter(BLOCK_DATA, index%in%indices)$hub_gene)
  # print(paste(nb_eQTL, nb_DE, nb_HUB))
}
names(TEST_RES)      = pheno_var
rownames(TEST_RES)   = Cluster
TEST_RES$Parallelism = Parallelism
TEST_RES$Parallelism = Parallelism/max(Parallelism)
TEST_RES             = TEST_RES[order(TEST_RES$Parallelism),]
corrplot(as.matrix(t(TEST_RES)))

############################
############################
############################

# indices
# 
# head(TEST_RES)
# 
# TO_TEST = c("clust_3", "evol_response_pool", "eQTL_carrier", "eQTL_phenotype", "is_ASE", "hub_gene", "significant_sg_pool")
# 
# ggcorrplot(RES$fisher_signif, type="upper")
# RES$fisher_signif["clust_3",-1]
# 
# fviz_dend(res)
# fviz_dist(res.dist)
# plot(res.hc, cex = 0.5)
# head(round(as.matrix(res.dist), 2))[, 1:20]
# sort(unique(as.vector(res.dist)))
# 
# max(cutree(res.hc, h=5))
# 
# res.hc$
# res.hc$height
# 
# KEYS = c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift")
# D = PERM_DATA[PERM_DATA$isShifting==1 & PERM_DATA$chr_int==1,]
# M = D[,KEYS]
# 
# fit = pvclust(t(M), method.hclust="ward", method.dist="binary")
# plot(fit) # dendogram with p values
# # add rectangles around groups highly supported by the data
# pvrect(fit, alpha=.99)
# table(cutree(fit$hclust, h=0.1))
# clusters <- pvpick(fit)
# clusters$clusters[[1]]
# 
# # as.numeric(clusters$clusters[[2]])
# D[D$index%in%c(54,58,41,42,11,15,2,62,22),KEYS]
# D$haplotype_block
# head(D)

############################
############################
############################

# MAP     = build_eQTL_associations(HB_DATA, PERM_DATA, eQTLs)
# 
# support_list = filter(MAP, Mx2_eQTL==1)$hb_eQTL
# run_permutation_test(MAP, "Mx2_pheno", "hb_eQTL", support_list, 1000)
# 
# support_list = filter(MAP, L1_pheno==1)$hb_eQTL
# run_permutation_test(PERM_DATA, "isShifting", "haplotype_block", support_list, 1000)
# 
# filter(MAP, isShifting_eQTL==1)
# mean(filter(PERM_DATA, eQTL_phenotype==1)$isShifting)
# mean(filter(MAP, L1_eQTL==1)$L1_pheno)
# 
# filter(MAP, L1_eQTL==1 & L1_pheno==1)
# filter(MAP, L2_eQTL==1 & L2_pheno==1)
# filter(MAP, L3_eQTL==1 & L3_pheno==1)
# filter(MAP, Mx1_eQTL==1 & Mx1_pheno==1)
# 
# filter(MAP, isShifting_eQTL==1 & isShifting_pheno==1)
# 
# g = graph.data.frame(MAP[,c("hb_eQTL", "hb_pheno")])
# plot(g, layout=layout_as_tree, vertex.size=10, label.cex=0.5)
# filter(MAP, hb_pheno==15)
# filter(MAP, hb_eQTL==2770)
# filter(MAP, hb_pheno==535)
# rownames(PERM_DATA) = PERM_DATA$haplotype_block
# pleio = sort(table(MAP$hb_eQTL))
# shift = PERM_DATA[as.character(names(pleio)),"isShifting"]
# 
# D = data.frame(pleio, shift)
# mean(filter(D, pleio>1)$shift)
# 
# filter(PERM_DATA, haplotype_block%in%names(pleio))$haplotype_block
# names(pleio)

############################
############################
############################

# mean(filter(PERM_DATA, eQTL_carrier==1)$evol_response_pool)
# mean(PERM_DATA$evol_response_pool)
# 
# mean(filter(PERM_DATA, eQTL_carrier==1)$isShifting)/mean(PERM_DATA$isShifting)
# mean(filter(PERM_DATA, eQTL_phenotype==1)$isShifting)/mean(PERM_DATA$isShifting)
# mean(filter(PERM_DATA, is_ASE==1)$isShifting)/mean(PERM_DATA$isShifting)
# mean(filter(PERM_DATA, hub_gene==1)$isShifting)/mean(PERM_DATA$isShifting)
# 
# mean(filter(PERM_DATA, eQTL_carrier==1)$evol_response_pool)/mean(PERM_DATA$evol_response_pool)
# mean(filter(PERM_DATA, is_ASE==1)$evol_response_pool)/mean(PERM_DATA$evol_response_pool)
# mean(filter(PERM_DATA, eQTL_phenotype==1)$evol_response_pool)/mean(PERM_DATA$evol_response_pool)
# mean(filter(PERM_DATA, hub_gene==1)$evol_response_pool)/mean(PERM_DATA$evol_response_pool)
# 
# dim(PERM_DATA)
# sum(PERM_DATA$is_ASE)
# sum(PERM_DATA$evol_response_pool)
# sum(PERM_DATA$eQTL_carrier)
# sum(PERM_DATA$eQTL_phenotype)
# sum(PERM_DATA$hub_gene)
# 
# TO_TEST = c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift")
# TO_TEST = c("mean_L1_AFC_abs", "mean_L2_AFC_abs", "mean_L3_AFC_abs", "mean_L5_AFC_abs", "mean_L6_AFC_abs", "mean_Mx1_AFC_abs", "mean_Mx2_AFC_abs")
# N       = length(TO_TEST)
# M_mean  = matrix(rep(0.0, N*N), ncol=N)
# for(chrom in seq(1,10))
# {
#   RES = run_matrix_tests(PERM_DATA, TO_TEST, chrom)
#   M_mean = M_mean+RES$spearman_signif
# }
# M_mean = M_mean/10
# ggcorrplot(M_mean, hc.order=T)
# 
# RES = run_matrix_tests(PERM_DATA, TO_TEST, 0)
# ggcorrplot(RES$pearson_signif, hc.order=T)

# 
# filter(PERM_DATA, eQTL_carrier==1)$Parallelism
# filter(PERM_DATA, eQTL_carrier==1)$isShifting
# 


#hist(PERM_DATA$Parallelism)
# names(HB_DATA)
# 
# GG_D = list("Shifting"=filter(PERM_DATA, isShifting==1)$haplotype_block,
#             "Partial"=filter(PERM_DATA, PartiallyParallel==1)$haplotype_block,
#             "Parallel"=filter(PERM_DATA, isParallel==1)$haplotype_block,
#             "DE"=filter(PERM_DATA, evol_response_pool==1)$haplotype_block)
#             #"eQTL"=filter(PERM_DATA, eQTL_carrier==1)$haplotype_block,
#             #"Pheno"=filter(PERM_DATA, eQTL_phenotype==1)$haplotype_block)
# ggvenn(GG_D)
# 
# mean(filter(PERM_DATA, Mx1_shift==1)$evol_response_pool)
# mean(PERM_DATA$evol_response_pool)
# 
# run_permutation_test(PERM_DATA, "evol_response_pool", filter(PERM_DATA, Mx1_shift==1)$index, 1000)
# 
# PERM_DATA = PERM_DATA[order(PERM_DATA$chr_int),]
# 
# par(mfrow=c(7,1), mar=c(0,4,0,0))
# plot(PERM_DATA$isShifting, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$isParallel, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$evol_response_pool, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$eQTL_carrier, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$eQTL_phenotype, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$hub_gene, type="h", col=PERM_DATA$chr_int)
# plot(PERM_DATA$significant_sg_pool, type="h", col=PERM_DATA$chr_int)
# 
# unique(PERM_DATA$chr_int)
# run_permutation_test(PERM_DATA, "hub_gene", filter(PERM_DATA, chr_int==10)$index, 1000)
# 
# unique(PERM_DATA$chr_int)
# tapply(PERM_DATA$Parallelism, PERM_DATA$chr, mean)
# 
# manhattan(PERM_DATA, chr="chr_int", bp="chr_index", p="Parallelism")
#names(PERM_DATA)
#D = PERM_DATA[PERM_DATA$chr_int==1,c("L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift")]
#D = PERM_DATA[,c("max_L1_AFC_abs", "max_L2_AFC_abs", "max_L3_AFC_abs", "max_L5_AFC_abs", "max_L6_AFC_abs", "max_Mx1_AFC_abs", "max_Mx2_AFC_abs")]

# PP = list()
# for(i in seq(1:10))
# {
#   D = PERM_DATA[PERM_DATA$chr_int==i,c("mean_L1_AFC_abs", "mean_L2_AFC_abs", "mean_L3_AFC_abs", "mean_L5_AFC_abs", "mean_L6_AFC_abs", "mean_Mx1_AFC_abs", "mean_Mx2_AFC_abs")]
#   M = cor(scale((D)))
#   P = cor.mtest(scale((D)))$p
#   p = ggcorrplot(M, hc.order=T, p.mat=P, sig.level = 0.05, insig = "blank", tl.cex=0)
#   PP[[i]] = p
# }
# do.call("grid.arrange", c(PP, ncol=4))
# 
# D = PERM_DATA[,c("mean_L1_AFC_abs", "mean_L2_AFC_abs", "mean_L3_AFC_abs", "mean_L5_AFC_abs", "mean_L6_AFC_abs", "mean_Mx1_AFC_abs", "mean_Mx2_AFC_abs")]
# M = cor(scale(D))
# P = cor.mtest(scale(D))$p
# ggcorrplot(M, hc.order=T, lab=T, p.mat=P, sig.level = 0.05, insig = "blank")
# 
# filter(PERM_DATA, eQTL_carrier==1 & L1_shift==1)
# filter(PERM_DATA, eQTL_carrier==1 & L2_shift==1)
# filter(PERM_DATA, eQTL_carrier==1 & L3_shift==1)
# filter(PERM_DATA, eQTL_carrier==1 & L5_shift==1)
# filter(PERM_DATA, eQTL_carrier==1 & L6_shift==1)

# filter(PERM_DATA, eQTL_carrier==1)$Parallelism
