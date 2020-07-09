# integration of CRISPRi data into model

# 14/09/2018
# expect to receive CRISPRi data from our partners towards the end of 2018
#
# I downloaded 3 datasets (2 CRISPR  + 1 shRNA) of gene-invalidation screens that are also for cancer cell lines, with or without drugs
# to see how CRISPR scores can be related to the results of simulations

# these are stored as subfolders in the folder "~/research/data/CRISPR/" 
#
# Phospho-ERK is a biomarker of response to a synthetic lethal drug combination of sorafenib and MEK inhibition in liver cancer_2018
# around 500 kinases and 10 essential genes with and without kinase (c-raf) inhibitor sorafenib
# 
# Scientific_Data_4_High-throughput RNAi screen for essential genes and drug synergistic combinations in colorectal cancer
# SW620 cell line (mad 26 others), primary (800 genes) and secondary (40) screens. 3 drugs, synergy scores
#
# Gene Essentiality Profiling Reveals Gene Networks and Synthetic Lethal Interactions with Oncogenic Ras_CELL2017
# 14 cell lines, 6 ras-mutant, full genome screen

#############
# The main issue before receiving data is to develop a pipeline to introduce gene inhibitions into the model
# Problem: many genes in screen, very few in model
# 
# Including gene-inhibitions in model, the different cases:
# 1) genes directly represented in the model: set to 0 initial value, 0 rates
# 2) genes that are intermediate nodes along edges (that represent indirect connections): cut the edge
# 3) genes internal to modules represented by 'abstract' nodes (such as homologous repair): reduce transition rate (by how much?) of abstract node
# 4) genes external to model: find path(s), calculate sign of effect, reduce transition rate scaled by distance (and # of paths?)

# currently best performing model is in "~/results/dnadam_atmatr_chk1mk2_hr_celldeathdeactivation"

# CRISPR list of genes
library(readr); library(plyr); library(dplyr); library(tidyr)
# load functions for finding paths etc
setwd("~/research/models/KRAS_DNA_repair_model/")
source("~/research/models/KRAS_DNA_repair_model/crispr_model_integration_functions.R")

crispri_gene_list_annotated <- read_tsv("~/research/data/CRISPR/crispr_target_list_pipeline/20180420-crispri-gene-list-annotated.tsv")
crispri_gene_list <- read_tsv("~/research/data/CRISPR/crispr_target_list_pipeline/20180420-crispri-gene-list.txt")
if (colnames(crispri_gene_list)!="HGNC") {
  crispri_gene_list <- rbind(colnames(crispri_gene_list), crispri_gene_list)
  colnames(crispri_gene_list) <- "HGNC"
}
# convert to array
crispri_gene_list <- crispri_gene_list$HGNC
# setdiff(crispri_gene_list_annotated$HGNC[!duplicated(crispri_gene_list_annotated$HGNC)],crispri_gene_list$HGNC)
# are 2 lists identical? yes

# load full network from Omnipath
pypath_directed_interactions <- read_tsv("~/research/models/pypath_omnipath/pypath_directed_interactions_HGNC_UNIPROT_databases_references.tsv")
# keep only signed interactions
if (sum(pypath_directed_interactions$interaction_directed_signed == 0)>0) {
  pypath_directed_interactions <- pypath_directed_interactions[!(pypath_directed_interactions$interaction_directed_signed == 0),]  
}

# load full network from SIGNOR
# SIGNOR_all <- read_tsv("~/research/models/pypath_omnipath/SIGNOR_all_data_20_02_18.tsv")
# SIGNOR_all <- SIGNOR_all[,c("ENTITYA", "EFFECT", "ENTITYB", "IDA", "IDB")]
# colnames(SIGNOR_all) <- c("source_hgnc", "interaction_directed_signed", "target_hgnc", "source_uniprot", "target_uniprot")
# # binary signs
# SIGNOR_all$interaction_directed_signed[grepl("down", SIGNOR_all$interaction_directed_signed)] <- as.numeric(-1)
# SIGNOR_all$interaction_directed_signed[grepl("up", SIGNOR_all$interaction_directed_signed)] <- as.numeric(1)
# # which are HGNC? (ie genes), have 2 columns to indicate
# HGNC_uniprot_fullgenome <- read_tsv("~/research/models/colon_cancer_map/NETWORK/HGNC_uniprot_fullgenome.tsv")
# SIGNOR_all$source_hgnc_status <- SIGNOR_all$source_hgnc; SIGNOR_all$source_hgnc_status[SIGNOR_all$source_hgnc %in% HGNC_uniprot_fullgenome$HGNC_gene_symbol] <- 1
# SIGNOR_all$source_hgnc_status[!SIGNOR_all$source_hgnc %in% HGNC_uniprot_fullgenome$HGNC_gene_symbol] <- 0
# SIGNOR_all$target_hgnc_status <- SIGNOR_all$target_hgnc; SIGNOR_all$target_hgnc_status[SIGNOR_all$target_hgnc %in% HGNC_uniprot_fullgenome$HGNC_gene_symbol] <- 1
# SIGNOR_all$target_hgnc_status[!SIGNOR_all$target_hgnc %in% HGNC_uniprot_fullgenome$HGNC_gene_symbol] <- 0
# write_tsv(SIGNOR_all, "~/research/models/pypath_omnipath/SIGNOR_all_simplified.tsv")
SIGNOR_all=read_tsv("~/research/models/pypath_omnipath/SIGNOR_all_simplified.tsv")

#################################
#################################
# There are 4 cases foreseeable concerning the relationship of genes in CRISPRi screen to the model
# 0) genes in the model
# 1) Genes (g) that are on the path connecting two nodes of the model (M): M_i -> g -> M_j (define max length of connecting paths)
# 2) Genes that are outside the model, affecting one or more nodes of the model: g -(...)-> M_i , M_i <- g -> M_j
# 3) Genes that are components of 'higher-level' nodes that represent multiple genes

# current best performing version of model
best_model <- read_tsv("~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/results/krasmodel_repair_one_variable/dnadam_atmatr_chk1mk2_hr_celldeathdeactivation/best_model.csv")
model_nodes <- data.frame(cbind(unique(as.vector(as.matrix(best_model[,c(1,3)]))),"")); colnames(model_nodes) <- c("model_name","HGNC")
# # HGNC names for model names...this is MANUAL, CHECK!!!
model_nodes$HGNC <- c(NA,"KRAS",NA,"CHEK1","MAPKAPK2",NA,NA,"CDC25B",NA,NA)
# # have model in HGNC form
# if(colnames(best_model[,1])!="source_hgnc" & colnames(best_model[,3])!="target_hgnc") {
#   best_model_hgnc <- best_model[,c(3,2,1)]; colnames(best_model_hgnc) <- c("source_hgnc","interaction_directed_signed","target_hgnc")
#   ########
#   # with "merge
#   best_model_hgnc$order <- as.numeric(rownames(best_model_hgnc))
#   best_model_hgnc <- merge(best_model_hgnc,model_nodes, by.x="source_hgnc", by.y="model_name",sort=F)
#   best_model_hgnc <- merge(best_model_hgnc,model_nodes, by.x="target_hgnc", by.y="model_name",sort=F)
#   best_model_hgnc$source_hgnc[!is.na(best_model_hgnc$HGNC.x)] <- best_model_hgnc$HGNC.x[!is.na(best_model_hgnc$HGNC.x)]
#   best_model_hgnc$target_hgnc[!is.na(best_model_hgnc$HGNC.y)] <- best_model_hgnc$HGNC.y[!is.na(best_model_hgnc$HGNC.y)]
#   best_model_hgnc <- best_model_hgnc[order(best_model_hgnc$order),c("source_hgnc","interaction_directed_signed","target_hgnc")]
#   rownames(best_model_hgnc) <- c()
# }
# write_csv(best_model_hgnc, "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/best_model_hgnc_SIF.csv")
best_model_hgnc <- read_csv("~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/results/krasmodel_repair_one_variable/best_model_hgnc_SIF.csv")

# unified column names, x0 being model nodes
colnames_table<-c("x4","i4","x3","i3","x2","i2","x1","i1","x0")

############
# Indirect connections BETWEEN MODEL nodes

# # all interactions in the model
# all_pairs <- expand.grid(model_nodes$HGNC[!sapply(model_nodes$HGNC, is.na)], model_nodes$HGNC[!sapply(model_nodes$HGNC, is.na)])
# # remove auto-loops?
# # all_pairs <- all_pairs[!(all_pairs$Var1 == all_pairs$Var2),]
# rownames(all_pairs)<-c(); 
# # convert dataframe columns from factor to chars (cf [] notation)
# all_pairs[] <- lapply(all_pairs, as.character); colnames(all_pairs) <- c("source_hgnc","target_hgnc")
# # extract whole paths and create table: source -> X -> target
# source_column_number <- 1; target_column_number<-3; df_large_network <- SIGNOR_all; df_small_network <- all_pairs
# source <- all_pairs$source_hgnc; target <- all_pairs$target_hgnc


##################################################
##################################################
###### Connections between MODEL nodes themselves

# all paths connecting model nodes with each other in one table, up to 4 step (3 intermediate) paths
source <- model_nodes$HGNC[!is.na(model_nodes$HGNC)]; target<-model_nodes$HGNC[!is.na(model_nodes$HGNC)]; 
source_column_number<-1; target_column_number<-3
# all paths connecting model nodes with each other in one table, up to 4 step (3 intermediate) paths
model_model_all_paths <- function_internal_paths_model_nodes(SIGNOR_all, source, target,
                                                             source_column_number,target_column_number, colnames_table)
if (!file.exists("~/research/data/CRISPR/model_model_all_paths.csv")) {
  write_csv(model_model_all_paths, "~/research/data/CRISPR/model_model_all_paths.csv")
}

##################################################
##################################################
###### Connections between SCREEN nodes and MODEL nodes

# all paths connecting SCREEN GENES to MODEL NODES in one table, up to 4 step (3 intermediate) paths
crispri_gene_list_except_model <- crispri_gene_list[!crispri_gene_list %in% model_nodes$HGNC[!is.na(model_nodes$HGNC)]]
# source <- crispri_gene_list_except_model; target<-model_nodes$HGNC[!is.na(model_nodes$HGNC)]
# screen_model_all_paths <- function_screen_model_paths(SIGNOR_all, source, target, source_column_number, target_column_number, colnames_table)
# if (!file.exists("~/research/data/CRISPR/screen_model_all_paths.csv")) {
#  write_csv(screen_model_all_paths, "~/research/data/CRISPR/screen_model_all_paths.csv")
# }

screen_model_all_paths=read_csv("~/research/data/CRISPR/screen_model_all_paths.csv")

#########################
#########################
# screen genes not amongs 1st, 2nd and 3rd neighbors
length(setdiff(crispri_gene_list_except_model, unique(function_collapse_df_column_to_vector(screen_model_all_paths[,c(1,3,5,7)]),
                                                      function_collapse_df_column_to_vector(model_model_all_paths[,c(3,5,7)]))))
# overlap
length(intersect(crispri_gene_list_except_model, unique(function_collapse_df_column_to_vector(screen_model_all_paths[,c(1,3,5,7)]),
                                                      function_collapse_df_column_to_vector(model_model_all_paths[,c(3,5,7)]))))

#########################
#########################
# check with igraph if my functions above work properly

# from my own table
source <- c("RAF1","KRAS"); target <-c("CDC25B","CDC25A")
own_paths <- function_screen_model_paths(SIGNOR_all, source, target, source_column_number, target_column_number, colnames_table)
own_paths <- own_paths[apply(apply(own_paths, 1, is.na),2,sum) < ncol(own_paths),]
if (sum(source %in% target) == 0) {
  own_paths <- own_paths[apply(apply(own_paths[,c(1,3,5,7,9)],1,duplicated,incomparables=NA),2,sum)==0,]
}
rownames(own_paths) <- c()

# igraph
library('igraph')
# LOAD entire network
# net <- as.matrix(SIGNOR_all[,c(1,3)]) # file containing the edge list in the form node1 node2 weight
graph <- graph.data.frame(as.matrix(SIGNOR_all[,c(1,3)]), directed=T, vertices=NULL) # E(graph)$weight=as.numeric(net[,3])
adj<-as.matrix(get.adjacency(graph)) # get.adjacency(graph,attr = "weight")
# GET connections
igr_paths <- all_shortest_paths(graph, source, target, mode="out"); igr_path_lists <- lapply(igr_paths$res, names)
# for same length paths: # igr_paths_df<-data.frame(do.call(rbind,lapply(igr_paths$res, names))) # do.call(rbind, l)
# if different lengths
igr_paths_df <- data.frame( t(sapply(igr_path_lists, "[", seq_len(max(sapply(igr_path_lists, length))))) )
igr_paths_df <- igr_paths_df[!duplicated(igr_paths_df),]; rownames(igr_paths_df) <- c()
# 
# igraph results identical with my own scripts
#########################

# some sample paths
SIGNOR_all=read_tsv("~/research/models/pypath_omnipath/SIGNOR_all_simplified.tsv")
colnames_table<-c("x4","i4","x3","i3","x2","i2","x1","i1","x0")
sample_path <- function_screen_model_paths(SIGNOR_all, c("SMAD4", "NCK1", "TRIP12", "CUBN", "AKR1B15"), 
                                           c("SMAD4","MAPKAPK2"), 1, 3, colnames_table)
colnames(sample_path) <- c("source","interaction","target")
a<-sample_path[,1:3]; b<-sample_path[,3:5];c<-sample_path[,5:7];d<-sample_path[,7:9]
colnames(a) <- c("source","interaction","target"); colnames(b) <- c("source","interaction","target"); 
colnames(c) <- c("source","interaction","target"); colnames(d) <- c("source","interaction","target"); 
sample_path_SIF <- bind_rows(a,b,c,d); sample_path_SIF<-sample_path_SIF[!duplicated(sample_path_SIF),];
sample_path_SIF <- sample_path_SIF[!colSums(apply(sample_path_SIF, 1, is.na))>0,]; 
sample_path_SIF<-sample_path_SIF[order(sample_path_SIF$target),]
# write_tsv(sample_path_SIF[,3:1], "../sample_path_SIF.tsv")


# summing effects
m<-t(apply(sample_path[,c(2,4,6,8)],1,as.numeric)); m[is.na(m)]<-1; net_inter <- apply(m, 1, prod)
net_inter=net_inter*( 1/(4-colSums(apply(sample_path[,c(2,4,6,8)], 1, is.na))) ) # [colSums(apply(sample_path, 1, is.na))>0]
sum_net_inter=sum(net_inter)/length(net_inter)

# sum net effects by last nodes
converging_net_effects <- as.data.frame(cbind(net_inter,match(sample_path$x1, unique(sample_path$x1)))) %>% group_by(V2) %>% 
  summarise_all(funs(sum))

#########################

# Matching simulation values with CRISPRi scores

# FC vs prolif rate
gamma1=1.05; gamma2=1; delta_gamma=gamma2-gamma1; t=10
x1_0=5; x2_0=95; r_0=x1_0/(x1_0+x2_0)
# FC
( (x1_0*2^(t*gamma1))/(x1_0*2^(t*gamma1) + x2_0*2^(t*gamma2)) )/r_0
# FC = 
1/(r_0 + (1-r_0)*2^(t*delta_gamma))
# r_0 and t from experiment; delta_gamma from simulations
# can plug in r_0 from data, delta_gamma from model simulations, 
# and look at correlations between FC_calc and FC_exp
# or calculate delta_gamma from experiment, from simulations and compare

# more directly: take x1_0 and x1_t, calculate exponent x1_t = x1_0*2^(gamma*t), 
# and compare it to simulation value of proliferation node, with the default value subtracted or divided by it
# gamma_exp*t = log2(x1_t/x1_0), compared to (prolif_perturb - prolif_default) or (prolif_perturb/prolif_default)


##############
# mapply examples
# mapply(function(x, y) seq_len(x) + y,
#        c(a=1, b=2, c=3),  # names from first
#        c(A=10, B=0, C=-10))
# 
# mapply(rep, "a", 1:4)