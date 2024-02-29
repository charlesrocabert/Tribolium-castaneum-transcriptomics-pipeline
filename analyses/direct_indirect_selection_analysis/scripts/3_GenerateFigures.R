#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_GenerateFigures.R
# -------------------
# Generate all figures associated to the manuscript.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("cowplot")
library("ggpubr")
library("ggbreak")
library("patchwork")
library("RColorBrewer")

parallel_category_figures <- function( gene_dataset, eQTL_dataset )
{
  ############
  comparisons = list(c("No significant shift", "Single shift"),
                     c("No significant shift", "Partial or full parallelism"),
                     c("Single shift", "Partial or full parallelism"))
  ############
  count_summary_1 = gene_dataset %>% group_by(Parallel_category) %>% tally()
  count_summary_2 = eQTL_dataset %>% group_by(Parallel_category) %>% tally()
  ############
  p1 = ggplot(gene_dataset, aes(Parallel_category, log_net_selection, fill=Parallel_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary_1$n), x=c(1,2,3), y=rep(3.5,3), size=3, vjust=2) +
    ylim(-6, 6.2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("Gene sequences") +
    labs(fill = "Parallelism:") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p2 = ggplot(eQTL_dataset, aes(Parallel_category, log_net_selection, fill=Parallel_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary_2$n), x=c(1,2,3), y=rep(1.6,3), size=3, vjust=2) +
    ylim(-6, 4) +
    xlab("") +
    ylab("Absolute weighted net selection\n(log-scale)") +
    ggtitle("eQTLs") +
    labs(fill = "Parallelism:") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p_no_legend = plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    ncol=2, labels="AUTO")
  legend = as_ggplot(get_legend(p1 + theme(legend.position="bottom", legend.text=element_text(size=10), legend.key.size = unit(1, 'cm'), legend.title=element_text(face="bold"))))
  ############
  p = plot_grid(p_no_legend, legend, ncol=1, rel_heights=c(1,0.1))
  ############
  return(p)
}

eQTL_figures <- function( eQTL_carrier_dataset, eQTL_phenotype_dataset )
{
  ############
  comparisons1 = list(c("Other", "eQTL carrier"))
  comparisons2 = list(c("Other", "eQTL phenotype"))
  # comparisons3 = list(c("No pleiotropy", "Low pleiotropy"),
  #                     c("No pleiotropy", "High pleiotropy"),
  #                     c("Low pleiotropy", "High pleiotropy"))
  comparisons3 = list(c("No pleiotropy", "Pleiotropy"))
  comparisons4 = list(c("Single eQTL", "Low connectivity"),
                      c("Single eQTL", "High connectivity"),
                      c("Low connectivity", "High connectivity"))
  # comparisons4 = list(c("Single eQTL", "At least two eQTLs"))
  ############
  count_summary_1 = gene_dataset %>% group_by(eQTL_carrier_category) %>% tally()
  count_summary_2 = gene_dataset %>% group_by(eQTL_phenotype_category) %>% tally()
  count_summary_3 = eQTL_carrier_dataset %>% group_by(Pleiotropy_category) %>% tally()
  count_summary_4 = eQTL_phenotype_dataset %>% group_by(Connectivity_category) %>% tally()
  ############
  p1 = ggplot(gene_dataset, aes(eQTL_carrier_category, log_net_selection, fill=eQTL_carrier_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons1, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary_1$n), x=c(1,2), y=rep(3,2), size=3, vjust=2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("eQTLs carriers") +
    labs(fill = "") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p2 = ggplot(gene_dataset, aes(eQTL_phenotype_category, log_net_selection, fill=eQTL_phenotype_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons2, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary_2$n), x=c(1,2), y=rep(3,2), size=3, vjust=2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("eQTLs phenotypes") +
    labs(fill = "") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p3 = ggplot(eQTL_carrier_dataset, aes(Pleiotropy_category, log_net_selection, fill=Pleiotropy_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons3, method="wilcox.test", label="p.signif") +
    #annotate("text", label=paste0("n = ",count_summary_3$n), x=c(1,2,3), y=rep(2.3,3), size=3, vjust=2) +
    annotate("text", label=paste0("n = ",count_summary_3$n), x=c(1,2), y=rep(2.3,2), size=3, vjust=2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("eQTL carriers pleiotropy") +
    labs(fill = "Pleiotropy:") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p4 = ggplot(eQTL_phenotype_dataset, aes(Connectivity_category, log_net_selection, fill=Connectivity_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons4, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary_4$n), x=c(1,2,3), y=rep(1.8,3), size=3, vjust=2) +
    #annotate("text", label=paste0("n = ",count_summary_4$n), x=c(1,2), y=rep(1.8,2), size=3, vjust=2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("eQTL phenotypes connectivity") +
    labs(fill = "Connectivity:") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  p = plot_grid(p1, p2, p3, p4, labels="AUTO")
  ############
  return(p)
}

hub_gene_figures <- function( gene_dataset )
{
  ############
  comparisons = list(c("Other","Hub gene"))
  ############
  count_summary = gene_dataset %>% group_by(hub_gene_category) %>% tally()
  ############
  p1 = ggplot(gene_dataset, aes(hub_gene_category, log_net_selection, fill=hub_gene_category)) +
    geom_boxplot() +
    scale_fill_brewer(palette="BrBG") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
    annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(3,2), size=3, vjust=2) +
    xlab("") +
    ylab("Absolute net selection\n(log-scale)") +
    ggtitle("Hub genes total selection compared to other genes") +
    labs(fill = "") +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  # p2 = ggplot(gene_dataset, aes(hub_gene_category, log10(Pleiotropy), fill=hub_gene_category)) +
  #   geom_boxplot() +
  #   scale_fill_brewer(palette="BrBG") +
  #   stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.signif") +
  #   annotate("text", label=paste0("n = ",count_summary$n), x=c(1,2), y=rep(2,2), size=3, vjust=2) +
  #   xlab("") +
  #   ylab("eQTL gene carrier pleiotropy\n(log-scale)") +
  #   ggtitle("Hub genes pleiotropy\ncompared to other genes") +
  #   labs(fill = "") +
  #   theme_classic() +
  #   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ############
  return(p1)
}

eQTL_histograms <- function( eQTL_carrier_dataset, eQTL_phenotype_dataset )
{
  ############
  p1 = ggplot(eQTL_carrier_dataset, aes(x=Pleiotropy)) +
    geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
    scale_y_cut(breaks=c(25), which=c(1), scales=c(0.5)) +
    xlab("Pleiotropy of eQTL carrier genes") +
    ylab("Count") +
    ggtitle("Distribution of eQTL carrier genes' pleiotropy") +
    theme_classic()
  ############
  p2 = ggplot(eQTL_phenotype_dataset, aes(x=Connectivity)) +
    geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
    scale_y_cut(breaks=c(100), which=c(1), scales=c(0.5)) +
    xlab("Phenotypes connectivity") +
    ylab("Count") +
    ggtitle("Distribution of phenotypes' connectivity") +
    theme_classic()
  ############
  #p1 = plot_grid(p1, labels=c("A"))
  #p2 = plot_grid(p2, labels=c("B"))
  #p  = plot_grid(p1, p2)
  return(list("p1"=p1, "p2"=p2))
  #return(p)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#------------------------------------------------#
# 1) Load the data                               #
#------------------------------------------------#
SNP_dataset            = readRDS("./analyses/direct_indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_dataset.rds")
eQTL_carrier_dataset   = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_carrier_dataset.rds")
eQTL_phenotype_dataset = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")

#------------------------------------------------#
# 2) Make parallelism figures                    #
#------------------------------------------------#
p = parallel_category_figures(gene_dataset, eQTL_dataset)
p
ggsave("analyses/direct_indirect_selection_analysis/plots/net_selection_parallelism.pdf", p, width=7, height=3.5, units="in")

#------------------------------------------------#
# 3) Make eQTL figures                           #
#------------------------------------------------#
# p = eQTL_figures(eQTL_carrier_dataset, eQTL_phenotype_dataset)
# ggsave("analyses/direct_indirect_selection_analysis/plots/net_selection_eQTLs.pdf", p, width=10, height=7, units="in")

#------------------------------------------------#
# 3) Make hub gene figures                       #
#------------------------------------------------#
# p = hub_gene_figures(gene_dataset)
# ggsave("analyses/direct_indirect_selection_analysis/plots/net_selection_hub_genes.pdf", p, width=5, height=4, units="in")

#------------------------------------------------#
# 4) Make pleiotropy and connectivity histograms #
#------------------------------------------------#
# figs = eQTL_histograms(eQTL_carrier_dataset, eQTL_phenotype_dataset)
# ggsave("analyses/direct_indirect_selection_analysis/plots/pleiotropy_distribution.pdf", figs[["p1"]], width=5, height=4, units="in")
# ggsave("analyses/direct_indirect_selection_analysis/plots/connectivity_distribution.pdf", figs[["p2"]], width=5, height=4, units="in")

