# this is a script to connect groups of nodes within a large network via indirect paths of up to 4 steps

library(readr); library(plyr); library(dplyr); library(tidyr)
# load SIGNOR network (PATH!!)
setwd("~/research/models/KRAS_DNA_repair_model/")
SIGNOR_all=read_tsv("network_paths/SIGNOR_all_simplified.tsv")
# load the functions to run the script, function file should be in same directory
source("network_paths/pathfinder_functions.R")
# names for columns
colnames_table=c("x4","i4","x3","i3","x2","i2","x1","i1","x0"); colnames_3cols=c("source","interaction","target")
source_column_number=1; target_column_number=3

##############################################################################
##############################################################################
# connecting nodes from Jane's regression analysis
# february/2019

# 2 groups of nodes we want to connect
source_nodes = c("SMAD4","NCK1","TRIP12","CUBN","AKR1B15"); 
target_nodes = c("KRAS","CHEK1","MAPKAPK2","CDC25B","ATM","ATR","FANCD2")

# RUN SCRIPT
sample_path <- function_screen_model_paths(SIGNOR_all, source_nodes, target_nodes, 1,3,colnames_table)
# colnames(sample_path) <- c("source","interaction","target")

# BRING it to SIF format
a<-sample_path[,1:3]; b<-sample_path[,3:5]; c<-sample_path[,5:7];d<-sample_path[,7:9]
colnames(a) <- colnames_3cols; colnames(b) <- colnames_3cols; colnames(c) <- colnames_3cols; colnames(d) <- colnames_3cols
sample_path_SIF <- bind_rows(a,b,c,d); sample_path_SIF<-sample_path_SIF[!duplicated(sample_path_SIF),];
sample_path_SIF <- sample_path_SIF[!colSums(apply(sample_path_SIF, 1, is.na))>0,]; 
sample_path_SIF<-sample_path_SIF[order(sample_path_SIF$target),]

# SAVE
# write_tsv(sample_path_SIF[,3:1], "network_paths/regression_to_model_paths_SIF.tsv")

# only paths of 2 steps max
sample_path_3steps = sample_path[,3:ncol(sample_path)]
# only where first column is a source node
sample_path_3steps=sample_path_3steps[sample_path_3steps$x3 %in% source_nodes,]
a<-sample_path_3steps[,1:3]; b<-sample_path_3steps[,3:5]; c<-sample_path_3steps[,5:7]
colnames(a)=colnames_3cols; colnames(b)=colnames_3cols; colnames(c) <- colnames_3cols

sample_path_SIF <- bind_rows(a,b,c); sample_path_SIF<-sample_path_SIF[!duplicated(sample_path_SIF),];
sample_path_SIF <- sample_path_SIF[!colSums(apply(sample_path_SIF, 1, is.na))>0,]; 
sample_path_SIF<-sample_path_SIF[order(sample_path_SIF$target),]
write_tsv(sample_path_SIF[,3:1], "network_paths/regression_to_model_paths_3steps_SIF.tsv")

##############################################################################
##############################################################################
# october/2019
# I will connect the nodes from Natalie's phosphoprot experiments and from the mutation data
# phosphoprot experiments data is contained in "~/KRAS_DNA_repair_model/data/wes_facs/wes_facs_data_tidy.csv" and its script is in 
# "KRAS_DNA_repair_model/phosphoprot_data.R"
# 
# mutation data is in 
# ~/KRAS_DNA_repair_model/data/mutation/mutation_data_tidy.csv
# ~/KRAS_DNA_repair_model/data/mutation/mutation_num_genes.csv
# and its script is also in "KRAS_DNA_repair_model/phosphoprot_data.R"

###########
# genes from mutation analysis
# source_nodes = c("SMAD4","NCK1","TRIP12","CUBN","AKR1B15");
mutation_num_genes = read_csv("data/mutation/mutation_num_genes.csv")
mutation_nodes=as.character(unique(mutation_num_genes$variable[mutation_num_genes$value>=4])[!grepl("KRAS|BRAF|MSI",unique(mutation_num_genes$variable[mutation_num_genes$value>=4]))])
mutation_nodes=gsub("\\(","",sapply(strsplit(mutation_nodes,"_"), "[[",1))
# all found in SIGNOR?
# length(intersect(c(unique(SIGNOR_all$source_hgnc),unique(SIGNOR_all$target_hgnc)), mutation_nodes))/length(mutation_nodes)
# yes

# model nodes
model_nodes = c("KRAS","CHEK1","MAPKAPK2","CDC25B","ATM","ATR","FANCD2")

###########
# nodes from mutation readouts
paths_mutations_model <- function_screen_model_paths(SIGNOR_all,mutation_nodes,model_nodes,1,3,colnames_table)
# paths with up to 3 intermediates
paths_mutations_model_SIF_three_intermediates = function_convert_SIF_three_intermediates(paths_mutations_model,colnames_3cols)
# with 2 intermediates: 87 edges
paths_mutations_model_SIF_two_intermediates=function_convert_SIF_two_intermediates(paths_mutations_model,mutation_nodes,colnames_3cols)
# with 1 intermediate: 14
paths_mutations_model_SIF_one_intermediate = function_convert_SIF_one_intermediate(paths_mutations_model,mutation_nodes,colnames_3cols)
# paths_mutations_model_SIF_two_intermediates[paths_mutations_model_SIF_two_intermediates$source %in% mutation_nodes,]

###########
# nodes from phospoprot readouts
# CDC2, CDC25C, CDK2, CHK2, P38
readout_nodes = c("CDK1", "CDC25C", "CDK2", "CHEK2", "MAPK14")
paths_readouts_model <- function_screen_model_paths(SIGNOR_all, readout_nodes, model_nodes, 1, 3, colnames_table)
# around 600 edges if we include 3-intermediate paths
paths_readouts_model_SIF_three_intermediates = function_convert_SIF_three_intermediates(paths_readouts_model,colnames_3cols)
# 2 intermediate paths: 195 edges
paths_readouts_model_SIF_two_intermediates=function_convert_SIF_two_intermediates(paths_readouts_model,readout_nodes,colnames_3cols)
# 1 intermediate: 27 edges
paths_readouts_model_SIF_one_intermediate=function_convert_SIF_one_intermediate(paths_readouts_model,readout_nodes,colnames_3cols)

###########
# paths between readout and mutation nodes
paths_readouts_mutations <- function_screen_model_paths(SIGNOR_all, readout_nodes, mutation_nodes, 1, 3, colnames_table)
# 3 intermediates: 1500 edges
paths_readouts_mutations_SIF_three_intermediates = function_convert_SIF_three_intermediates(paths_readouts_mutations,colnames_3cols)
# 2 intermediates: 431 edges
paths_mutations_readouts_SIF_two_intermediates=function_convert_SIF_two_intermediates(paths_readouts_mutations,readout_nodes,colnames_3cols)
# 1 intermediate: 72 edges
paths_mutations_readouts_SIF_one_intermediate=function_convert_SIF_one_intermediate(paths_readouts_mutations,readout_nodes,colnames_3cols)

# concatenate the 1-intermediate paths
one_intermediate_paths=bind_rows(paths_mutations_model_SIF_one_intermediate,
                                 paths_readouts_model_SIF_one_intermediate,
                                 paths_mutations_readouts_SIF_one_intermediate)
one_intermediate_paths=one_intermediate_paths[order(one_intermediate_paths$source),]
write_csv(one_intermediate_paths,"network_paths/one_intermediate_paths.csv")

# SIF file of model
krasmodel15vars_SIF=read_csv("network_paths/krasmodel15vars.csv")

#########################################
# 29/October/2019: I received an influence graph from Natalie, built by Markus
# I implemented this in the cytoscape file "model_expansion.cys" (same folder)
# and as SIF: 
#
# MaBoSS simulations of the logical model "markus_network_edited_modeling_with_autoloops.bnet"
# showed CDC25B is never activated, so it would need to be connected to PI3K or MAPK elements so that it has a driver

charitemodel=read_csv("network_paths/markus_network_edited_modeling_with_autoloops.sif")
charitemodel_nodes=unique(function_collapse_df_column_to_vector(charitemodel[,c(1,3)]))

source_nodes=c("KRAS","NRAS","HRAS","BRAF","MAPK3","MAPK1","MAPK11","MAPK12","MAPK13","MAPK14","MAPKAPK2","JNK1","JNK2","JNK3","MAP2K2","PI3K","AKT")
# PATHS to CDC25B
paths_to_cdc25b <- function_screen_model_paths(SIGNOR_all, source_nodes,"CDC25B",1,3,colnames_table)
# paths_readouts_model_SIF_one_intermediate=function_convert_SIF_one_intermediate(paths_to_cdc25b,"CDC25B",colnames_3cols)
# paths to CDC25A
paths_to_cdc25a <- function_screen_model_paths(SIGNOR_all, source_nodes,"CDC25A",1,3,colnames_table)
# paths to CDC25C
paths_to_cdc25c <- function_screen_model_paths(SIGNOR_all, source_nodes,"CDC25C",1,3,colnames_table)

# MK2 inputs
paths_to_mk2 <- function_screen_model_paths(SIGNOR_all, source_nodes[!source_nodes %in% "MAPKAPK2"],"MAPKAPK2",1,3,colnames_table)