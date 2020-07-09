# MaBoss_wrapper_plotter_R_scripts
# (author: mihaly.koltai@curie.fr)

# below are functions to 
# introduce mutations/perturbations in MaBoSS models easily
# define initial states from simulations
# run simulations of multiple models in one go, create separate folders for them
# create plots from multiple MaBoss simulations and save them in their respective folders

# commenting and examples in process...

# prerequirements: 
# install MaBoSS 2.0 from http://maboss.curie.fr
# source Maboss environment: "source MaBoSS.env" on command line
# load all functions in the function file (some of them nested)
setwd("~/research/models/KRAS_DNA_repair_model/maboss_models")
source("../maboss_wrapper_plotting_functions.R")

# BOOLNET -> SIF
# bnet_table=function_extract_bnet_convert_to_sif("~/research/models/KRAS_DNA_repair_model/matlab_ode/maboss_analytical/model_files/breast_cancer_zanudo2017.bnet")
# write_csv(bnet_table,"./bnet_sif_table.csv")
# function_extract_bnet_convert_to_sif("/data/users/mkoltai/research/models/KRAS_DNA_repair_model/network_paths/markus_network_edited_modeling_with_autoloops.bnet")
# path_bnet="~/research/models/KRAS_DNA_repair_model/matlab_ode/maboss_analytical/model_files/"; name_bnet="krasmodel15vars.bnet"
# sif_table=function_extract_bnet_convert_to_sif(paste(path_bnet,name_bnet,sep =""))
# path_sif=paste(path_bnet,"other/",sep=""); name_sif=gsub('bnet','csv',name_bnet)
# write_csv(sif_table,paste(path_sif,name_sif,sep =""))

##########################################################
# Ensemble modeling: create all models by defining possible rules per node

# CREATE MODEL TOPOLOGIES
# run this script from from "~/research/models/KRAS_DNA_repair_model/maboss_models"

# list of variables from SIF file
SIF_input <- merged_kras_dna_dam <- read_tsv("~/research/models/KRAS_DNA_repair_model/first_other_models/dna_repair/merged_kras_dna_dam.sif")
variable_names=unique(as.vector(as.matrix(SIF_input[,c(1,3)])))

# create inputs from SIF table
list_inputs <- lapply(unique(SIF_input$target_hgnc), 
                      function(x) {
                        inter_signs_inhib_negation <- SIF_input$interaction_directed_signed[SIF_input$target_hgnc %in% x]==(-1)
                        inputs<-SIF_input$source_hgnc[SIF_input$target_hgnc %in% x]
                        inputs[inter_signs_inhib_negation] <- paste("!",inputs[inter_signs_inhib_negation],sep="")
                        paste(inputs, collapse = " ")
                        })
names(list_inputs) <- paste(unique(SIF_input$target_hgnc),"_inputs",sep="")

######################

# logical rules (currently manually defined)
list_all_rules <- list(
     cc_inputs="cc",
     KRAS_inputs="KRAS",
     DSB_inputs="(DSB | KRAS) & !(FAHRR | HRR2 | NHEJ)",
     CHEK1_inputs=c("ATM | ATR", "ATM & ATR"),
     MAPKAPK2_inputs=c("(ATM | ATR) & KRAS", "ATM & ATR | KRAS", "ATM & ATR & KRAS"),
     FAHRR_inputs="(FAHRR | DSB | FANCD2I) & !NHEJ & !cell_death",
     HRR2_inputs="(HRR2 | DSB | NHEJ) & !NHEJ & !cell_death & !FAHRR",
     CDC25B_inputs="(cc | KRAS) & (!CHEK1 & !MAPKAPK2) & !cell_death",
     g2m_trans_inputs="g2m_trans | CDC25B",
     cell_death_inputs="cell_death | (DSB & g2m_trans)",
     ATM_inputs=c("ATR | FAcore | DSB", "DSB", "(ATR | FAcore) & DSB"),
     ATR_inputs="ATM",
     FAcore_inputs="ATM | ATR",
     FANCD2I_inputs="ATM | ATR | DSB | FAcore",
     NHEJ_inputs="(DSB | NHEJ) & !FAHRR & !HRR2 & !cell_death",
     TP53_inputs="ATM | ATR | NHEJ")

# check if varnames and list of rules match
if (sum(variable_names %in% gsub("_inputs","",names(list_all_rules)))/length(variable_names)!=1) {
  print("variables and rules don't match!!")
} else { print("OK, variables and rules match!!") }

# create a directory where results will be eventually moved
results_dirname <- "krasmodel_repair_one_variable/best_model"
if (!dir.exists(paste("kras_model_topologies/results/",results_dirname,sep=""))) {
  dir.create(paste("kras_model_topologies/results/",results_dirname,sep=""))
}

# library(combinat)
# GENERATE all possible models (without looping)
all_models_inputs <- do.call(expand.grid, list_all_rules); colnames(all_models_inputs) <- gsub("_inputs","",names(list_all_rules))
# align var list with colnames of all_models_inputs
if (sum(variable_names == colnames(all_models_inputs))/ncol(all_models_inputs)!=1) { variable_names <- colnames(all_models_inputs) }
# with indices of rules
all_models_inputs_numeric <- do.call(expand.grid,lapply(lapply(lapply(list_all_rules,as.factor), as.numeric),sort));
colnames(all_models_inputs_numeric) <- gsub("_inputs","",names(list_all_rules))

metafolder <- "/bioinfo/users/mkoltai/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies"
subfolder_stem <- 'krasmodel'

# save all models in a text file
write_csv(all_models_inputs, paste(metafolder,"all_models_inputs.csv",sep="/"))
write_csv(all_models_inputs_numeric, paste(metafolder,"all_models_inputs_numeric.csv",sep="/"))

#####################
# GENERATE all bnet, bnd and config (cfg) FILES + shell scripts to run them

# PARAMETERS for config files
# names for cell types (here) defined by initial conditions
filename_tag_list <- c("_wt","_DNAdam","_KRASmut")
# initstate lists: the nodes whose initial states are defined
nodes_to_set_initstate <- variable_names[variable_names %in% c("cc","KRAS","DSB")]
# the values these nodes take on in different versions of model
initstate_values_list <- list(c(1,0,0), c(1,0,1), c(1,1,1)); names(initstate_values_list)=c("cc","KRAS","DSB")
# FIXED GENES
# to what values these genes would be fixed: -1 is KO, +1 knock-in
fixed_genes <- variable_names[variable_names %in% c("CHEK1","MAPKAPK2")]
indices_fixed_genes <- list(c(-1,0),c(0,-1),c(-1,-1))

#####################

visible_nodes <- variable_names[!variable_names %in% c("cc","phosphatase","TP53")]
meta_sh_file_path <- "kras_model_topologies/run_all.sh"; path_bnet_to_maboss_py <- ""
maboss_path <- "~/research/models/MaBoSS-env-2.0/"

# this generates the mutation+init cond specific cfg files in all the model-specific subfolders, 
# and a sh script in all those subfolders
# "bnet_to_maboss_py" sets max time, threads, # of trajectories
all_models_inputs_matrix <- t(all_models_inputs)
function_bnet_bnd_cfg_model_creator(all_models_inputs_matrix, 
                      paste(gsub(paste(getwd(),"/",sep=""),"",metafolder),subfolder_stem,sep="/"), path_bnet_to_maboss_py, 
                      variable_names, visible_nodes,
                      filename_tag_list, nodes_to_set_initstate, initstate_values_list, 
                      fixed_genes, indices_fixed_genes,maboss_path, meta_sh_file_path)

# all cfg files
input_cfg_files <- system("ls kras_model_topologies/krasmodel*/*[a-z].cfg",intern=T)
#
# change UP or DOWN transition rates
# DNA damage repair slower than phosphoryl events
if(sum(variable_names %in% "DSB")>0) { function_change_rates(input_cfg_files, variable_names[variable_names %in% "DSB"], NA, 0.15) }
# CDC25B activation slower
if(sum(variable_names %in% "CDC25B")) { function_change_rates(input_cfg_files, variable_names[variable_names %in% "CDC25B"], 0.25, 1) }
# repair processes should be shut down immediately if cell death happens
lapply(c("FAHRR", "HRR2", "NHEJ"), function(x) { # "FAHRR", "HRR2", "NHEJ"
if(sum(variable_names %in% x)) { function_change_rates(input_cfg_files, x, NA, 100) } } )

# rewrite meta sh file if it already existed
subfolder_dirs <- function_subfolder_string(metafolder,subfolder_stem)
meta_sh_string <- function_meta_sh_string(subfolder_stem,subfolder_dirs,"/run_simul.sh")
write(meta_sh_string, "kras_model_topologies/run_all.sh")

# RUN ALL from command line (go to "kras_model_topologies"): 
# sh run_all.sh
#
# remove sub-subfolders: 
# rm -r -d krasmodel*/*/

##########################################################
##########################################################

# EXTRACT SIMULATION RESULTS

# run functions from "~/research/models/KRAS_DNA_repair_model/maboss_models"

######
# create lists of folder names
#
max_time=29.98; timevector <- seq(0,max_time,0.01); timestep_num=length(timevector)

all_cond_subfolders <- paste("/",list.dirs(subfolder_dirs,recursive=F),sep="")
# desired order of conditions (as they are automatically in all_cond_subfolders)
conditions_order <- c(5,6,8,7, 9,10,12,11, 1,2,4,3)

# reorder&rename subfolders (careful: hardwired variables!)
all_cond_subfolders <- function_reorder_allcond_folders(subfolder_dirs, all_cond_subfolders, conditions_order, timestep_num)

# the order of the conditions (cell lines + drugs) how they are shown on plots. This is very important otherwise all plots mixed up
# in later commands (plotting) this order should not be changed
conditions_list <- unique(gsub("model[0-9]+_","", unlist(lapply(
  sapply(strsplit(all_cond_subfolders,"/"),`[`,unique(sapply(strsplit(all_cond_subfolders,"/"), length))),
  function(x) rep(x,timestep_num)))))
conditions_list <- gsub("chek1off","CHEK1off",gsub("mk2off","MAPKAPK2off",conditions_list))
# conditions_list = c("wt","wt_fixed_chek1off","wt_fixed_mk2off","wt_fixed_chek1off_mk2off",
# "DNAdam","DNAdam_fixed_chek1off","DNAdam_fixed_mk2off","DNAdam_fixed_chek1off_mk2off",
# "KRASmut","KRASmut_fixed_chek1off","KRASmut_fixed_mk2off","KRASmut_fixed_chek1off_mk2off")

######
# CREATE (empty) DF to store results in
# hardwired: column names "model_folder", "condition_folder"
# visible_nodes <- variable_names[!variable_names %in% c("cc","phosphatase","TP53")]
column_names<-c("Time", "model_folder", "condition_folder", visible_nodes)
all_simulation_timecourses_df <- function_create_simul_empty_df(all_cond_subfolders, rownumber,timestep_num, 
                                                                timevector, column_names, "krasmodel", "model[0-9]+_")

######
# these parameters needed for internal functions of extraction
font_size <- 7; plot_width <- 6; plot_height <- 4; font_size<-12; plotting_threshold <- 0.001; line_thickness <- 1
# EXTRACT DATA of all simulations (round to 3 decimals)
ptm <- proc.time()
directory_depth <- 1:3; directory_index<-2
all_simulation_timecourses_df <- function_extract_simuls(all_cond_subfolders, all_simulation_timecourses_df, visible_nodes, 
                                                         timestep_num, directory_depth, directory_index)
proc.time() - ptm

# save
write_csv(all_simulation_timecourses_df,"kras_model_topologies/all_simulation_timecourses_df.csv")

# memory usage in MByte
# sizes_megabyte <- sort(unlist(sapply(ls(), function(n) object.size(get(n)), simplify=FALSE)),decreasing=T)/2^20
# sizes_megabyte[sizes_megabyte>0.5]
########################

# or READ IN FROM EXISTING FILE
# loaded_results_filepath <- "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/results/krasmodel_repair_one_variable/dnadam_chk1mk2_atmatr/all_simulation_timecourses_df.csv"
# all_simulation_timecourses_df <- read_csv(loaded_results_filepath)

# MAP ACTIVITY to PHOSPHORYLATION
# Boolean simulations describe activity, but we want to fit it to PHOSPHORYLATION levels, which is sometimes the inverse
vars_to_invert_activ_phos <- "CDC25B"
all_simulation_timecourses_df[is.na(all_simulation_timecourses_df)] <- 0
all_simulation_timecourses_df[,vars_to_invert_activ_phos] <- 1 - all_simulation_timecourses_df[,vars_to_invert_activ_phos]

# endpoints
all_simulation_endpoints <- all_simulation_timecourses_df[all_simulation_timecourses_df$Time %in% 
                                                      max(all_simulation_timecourses_df$Time),]
# all_simulation_endpoints[is.na(all_simulation_endpoints)]<-0

######
# max points for SOME variables
max_variables=c("CHEK1","MAPKAPK2")
all_simulation_maxpoints <- function_maxpoints(all_simulation_timecourses_df, all_simulation_endpoints, visible_nodes,max_variables, "Time")
# all_simulation_maxpoints$variable <- as.factor(all_simulation_maxpoints$variable)
# visible_nodes <- colnames(all_simulation_endpoints)[which(colnames(all_simulation_endpoints) %in% "kras"):ncol(all_simulation_endpoints)]

##########################################################
# FITTING/MODEL SELECTION - euclidean distance of simuls and data

# CAREFUL: "variable_names" or "visible_nodes"

# read in and structure data as tidy data
exp_values <- function_format_expvalues(as.data.frame(read_csv("exp_values_heatmap_rownames.csv",col_names = F)),"condition_folder",conditions_list)
exp_values$condition_folder=gsub("chek1off","CHEK1off",gsub("mk2off","MAPKAPK2off",exp_values$condition_folder))
# or read in from saved tidy df
# exp_values <- read_csv("~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/exp_values_tidy.csv")

# CALCULATE MSE
mse <- function_calc_mse(all_simulation_maxpoints, exp_values)

#################################################
# SAVE best performing models
# all_models_inputs[order(mse)[1:smallest_n],]
smallest_n <- 4
name_best5_sif <- "best5_models_redundant.sif"; 
metafolder <- "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies"
function_save_best_models(all_models_inputs,smallest_n,mse,name_best5_sif,metafolder)

##########################################################
##########################################################
# PLOTTING

# PLOT MSE values, smallest n values highlighted
smallest_n <- min(c(3,length(mse)))
postscript("kras_model_topologies/bestmodel_indices_mse.eps",height=height_width_vals[1],width=height_width_vals[2], 
           onefile=FALSE, horizontal=FALSE, paper='special')
n_to_show <- min(c(smallest_n*10, length(mse))); title_mse="MSE"; title_mse="model index"
function_mse_plot(mse, n_to_show,smallest_n,title_mse,xlab_var)
dev.off()
##################################

# PLOT DYNAMICS for 1 model topology, subplot with 3 rows: celltypes, 4 columns: conditions
# 1 model, 1 or more conditions

model_id <- order(mse)[1] # condition_id <- 1:length(conditions_list)
# extract data points
vars_to_plot <- c("MAPKAPK2","CHEK1","CDC25B","DSB","cell_death","g2m_trans") # "cell_death","DSB","CDC25B","MAPKAPK2", 
  # c("cell_death","DSB","g2m_trans",variable_names[grepl("HRR",variable_names)],"NHEJ","CDC25B","MAPKAPK2", "CHEK1") 
  # visible_nodes

time_limit <- 20
cond_folder_colname="condition_folder"
indiv_model_dynamics <- function_indiv_model_dynamics(all_simulation_timecourses_df,model_id,"model_folder",cond_folder_colname,
                                                      conditions_list,vars_to_plot,"Time")

# PLOT DYNAMICS under diff. conditions as separate subplots ("facets")
height_width_vals <- c(6,12);
postscript(paste("kras_model_topologies/dynamics_model",model_id,".eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2])
#######
ggplot(indiv_model_dynamics, aes(x=Time,y=value,color=variable),xlab="") + geom_line(size=1.5) +
  facet_wrap(~factor(condition_folder,levels=conditions_list),ncol=4) + # [match(1:length(conditions_list),conditions_order)]
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title=element_blank(),
                     legend.box.background=element_rect(), legend.box.margin=margin(1,1,1,1) )
#######
dev.off()

#######################################################
# PLOT MAXPOINTS of models without data being shown (show error or values)

title_indices <- as.character(1:length(unique(all_simulation_maxpoints$model_folder)))
# LOAD DATA, which VARS
vars_to_plot <- c("CHEK1", "MAPKAPK2", "DSB", "CDC25B", "cell_death") # , "proliferation"
model_id <- order(mse)[1] # order(mse)[1] unique(all_simulation_maxpoints$model_folder)
# select those model versions we want to plot
df_plot <- all_simulation_maxpoints[all_simulation_maxpoints$variable %in% vars_to_plot # variable to plot
                                    & all_simulation_maxpoints$model_folder %in% model_id,] # models to plot
df_plot$mse <- paste("MSE=", round(mse[df_plot$model_folder],3),sep="")
if (class(df_plot$value)!="numeric") { df_plot$value <- as.numeric(df_plot$value) }
df_plot <- merge(df_plot, exp_values, by = c("variable","condition_folder"))

# PLOT PARAMETERS
my_lines <- data.frame(x=c(4.5,8.5), y=c(0.5,0.5), xend=c(4.5,8.5), yend=rep(length(vars_to_plot)+0.5,2) )
# theme settings
x_y_size <- 10
strip_text_x_size <- x_y_size; axis_text_y_size <- x_y_size; legend_text_size <- 9; geom_text_size<-4; height_width_vals <- c(7,12); 
# theme1 defined in fcn file
# how many columns
col_num <- round(sqrt(length(model_id)))
# CHOOSE: plotting error or value!!!
what_to_plot_list <- c("vals","round_vals","errs","round_errs"); what_to_plot<-what_to_plot_list[1]
df_plot <- function_col_select(df_plot, what_to_plot)

# SAVE
postscript(paste("kras_model_topologies/heatmaps_","model_", paste(as.character(model_id),collapse = "_"),what_to_plot,".eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],
           onefile=FALSE, horizontal=FALSE, paper='special')
#######
# launch GGPLOT
ggplot(df_plot, aes(x=factor(condition_folder,levels=conditions_list), # conditions_list[match(1:length(conditions_list),conditions_order)]
      y=factor(variable, levels=rev(variable_names)), fill=value )) + # abs_dev_data value.x
  scale_fill_gradient(low="white",high="red",limits=c(0,1)) +
  geom_tile(color="black",size=0.1) + facet_wrap(~model_folder+mse,ncol=col_num,labeller=label_wrap_gen(multi_line=FALSE)) +
  geom_text(aes(label=round(value,1)), size=geom_text_size) + # values on plot
  theme1 + # theme settings
  geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.25,inherit.aes=F,color="blue")
#######
dev.off()

#######################################################
# PLOT MAXPOINT of selected models NEXT to DATA

theme1$axis.text.y$size<-13; theme1$strip.text.x<-element_text(size=15); height_width_vals <- c(5,12); 
theme1$panel.spacing=unit(2,"lines"); theme1$strip.background=element_rect(colour="white",fill="white"); theme1$panel.border=element_blank()
#######
# model to show, number of col.s
model_id <- order(mse)[1]; col_num <- 2
# VAR.s to plot
vars_to_plot <- c(unique(exp_values$variable)) # 'proliferation'
# LOAD DATA
# if (class(all_simulation_maxpoints$value)!="numeric") { all_simulation_maxpoints$value <- as.numeric(all_simulation_maxpoints$value) }
df_plot <- function_df_indiv_heatmap(all_simulation_maxpoints,vars_to_plot,model_id,exp_values)
what_to_plot_list <- c("vals","round_vals")
####
# WHAT TO PLOT: vals or rounded vals
what_to_plot<-what_to_plot_list[2]
####
if (grepl("round",what_to_plot)) {df_plot$value <- round(df_plot$value)}

# create eps file
postscript(paste("kras_model_topologies/heatmap_model",paste(model_id,collapse="_"),what_to_plot,"_with_expval",".eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
#######
ggplot(df_plot, aes(x=factor(condition_folder,levels=conditions_list), # [match(1:length(conditions_list),conditions_order)]
      y=factor(variable, levels=rev(variable_names)), fill=value)) + 
  facet_wrap(~mse+model_folder, ncol=col_num, labeller=label_wrap_gen(multi_line=FALSE) ) +
  scale_fill_gradient(low="white",high="red", limits=c(0,1), na.value="transparent") + geom_tile(color="black",size=0.1) +
  geom_text(aes(label=round(value,1)), size=4) + theme1 + geom_segment(data=my_lines,aes(x,y,xend=xend,yend=yend),size=1.25,inherit.aes=F,color="blue")
#######
dev.off()

# best model SIF
# move files
if (dir.exists(paste("kras_model_topologies/results/",results_dirname,sep=""))){
  files_to_copy <- setdiff(list.files("kras_model_topologies",include.dirs = FALSE,recursive = F,full.names = T),
                           list.dirs("kras_model_topologies",recursive = F))
  files_to_copy <- files_to_copy[!grepl("exp_values.csv",files_to_copy)]
  
  file.copy(from=files_to_copy, to=paste("kras_model_topologies/results/",results_dirname,sep=""), 
            overwrite = FALSE, recursive = FALSE, copy.mode = TRUE)
  file.remove(files_to_copy)
  
}

############################################################
############################################################
# # adding DNA repair components
# # the DNA repair process has to be more detailed
# # I use the model https://www.ncbi.nlm.nih.gov/pubmed/26385365
# # to include different types of DNA repair
# 
# # this is currently the best model from the optimization above, with HGNC names
# best_model_hgnc <- read_csv("~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/results/krasmodel_repair_one_variable/best_model_hgnc_SIF.csv")
# 
# # Rodriguez model in BoolNet
# rodriguez_2015_bnet <- read_csv("~/research/models/KRAS_DNA_repair_model/dna_repair/rodriguez_2015.bnet")
# rodriguez_2015_sif <- function_boolnet_to_SIF(rodriguez_2015_bnet)
# # save
# write_tsv(rodriguez_2015_sif, "~/research/models/KRAS_DNA_repair_model/dna_repair/rodriguez_2015.sif")
# 
# # Rodriguez model starts with ICL leading to ADD and DSB.
# # I removed CHKREC, ICL, ADD, NUC1, NUC2, as these are processes that precede DSB, but we want to start from the condition of DSB. 
# # this should be discussed if DSB is a good starting point
# nodes_remove <- c("CHKREC","ICL","ADD","NUC1","NUC2","TLS")
# rodriguez_2015_sif_reduced <- rodriguez_2015_sif[!(rodriguez_2015_sif$factors %in% nodes_remove | rodriguez_2015_sif$targets %in% nodes_remove),]
# rodriguez_2015_sif_reduced <- rodriguez_2015_sif_reduced[order(rodriguez_2015_sif_reduced$targets, rodriguez_2015_sif_reduced$factors),]; 
# rownames(rodriguez_2015_sif_reduced) <- c()
# # length( unique(function_collapse_df_column_to_vector(rodriguez_2015_sif_reduced[,c(1,3)])) )
# 
# # index those variables present in current kras model (best_model_hgnc)
# rodriguez_2015_sif_reduced$origin <- "rodriguez"
# # nodes present or not in KRAS model
# rodriguez_2015_sif_reduced$kras_model_presence <- 0
# # these nodes are represented in KRAS model, although merged into 'meta-nodes'
# nodes_KRAS_model <- c("ATM","ATR", "FANCD2I", "HRR2", "FAHRR", "DSB")
# rodriguez_2015_sif_reduced$kras_model_presence[rodriguez_2015_sif_reduced$factors %in% nodes_KRAS_model |
#                                                rodriguez_2015_sif_reduced$targets %in% nodes_KRAS_model] <- 1
# # use HGNC name
# # use HGNC name
# rodriguez_2015_sif_reduced[rodriguez_2015_sif_reduced=="P53"] <- "TP53"
# colnames(rodriguez_2015_sif_reduced)[1:ncol(best_model_hgnc)] <- colnames(best_model_hgnc)
# # list of nodes
# nodes_rodriguez_2015_sif_reduced <- unique(function_collapse_df_column_to_vector(rodriguez_2015_sif_reduced[,c(1,3)]))
# 
# # now this model needs to be merged with kras model
# # use the same names
# best_model_hgnc_tomerge <- best_model_hgnc
# best_model_hgnc_tomerge[best_model_hgnc_tomerge=="dna_dam"] <- "DSB"
# best_model_hgnc_tomerge[best_model_hgnc_tomerge=="atm_atr"] <- "ATM,ATR"
# best_model_hgnc_tomerge[best_model_hgnc_tomerge=="hr"] <- paste(unique(function_collapse_df_column_to_vector(rodriguez_2015_sif_reduced[,c(1,3)]))[
#               grepl("HRR",unique(function_collapse_df_column_to_vector(rodriguez_2015_sif_reduced[,c(1,3)])))],collapse=",")
# best_model_hgnc_tomerge$origin <- "krasmodel"
# # split by commas for doubled nodes
# z <- best_model_hgnc_tomerge %>% mutate(source_hgnc = strsplit(source_hgnc, ",")) %>% unnest(source_hgnc) %>% 
#   mutate(target_hgnc = strsplit(target_hgnc, ",")) %>% unnest(target_hgnc)
# best_model_hgnc_tomerge <- z[,colnames(best_model_hgnc_tomerge)]
# # remove interactions between nodes of DNA model, except self-activation that stand for states that keep their state
# best_model_hgnc_tomerge <- best_model_hgnc_tomerge[
#                           !colSums(apply(best_model_hgnc_tomerge, 1, function(x) x %in% nodes_rodriguez_2015_sif_reduced ))==2 | # nodes from DNA model
#                           colSums(apply(best_model_hgnc_tomerge,1,duplicated))>0,] # nodes that have duplicates
# 
# merged_kras_dna_dam <- rbind(best_model_hgnc_tomerge, rodriguez_2015_sif_reduced[,1:ncol(best_model_hgnc_tomerge)])
# merged_kras_dna_dam <- merged_kras_dna_dam[!duplicated(merged_kras_dna_dam[,1:3]),]
# 
# write_tsv(merged_kras_dna_dam[,c(3,2,1,4)],"~/research/models/KRAS_DNA_repair_model/dna_repair/merged_kras_dna_dam.sif")
# # merged_kras_dna_dam=read_tsv("~/research/models/KRAS_DNA_repair_model/dna_repair/merged_kras_dna_dam.sif")

############################################################
############################################################

# PARAMETER SAMPLING: with modified MaBoSS scripts

# INPUT for parameter sampling w multiple conditions: generate JSON file
# CONDITIONS per cell line
## names for cell types (here) defined by initial conditions
# filename_tag_list <- c("_wt","_DNAdam","_KRASmut")
# initstate lists: the nodes whose initial states are defined, nodes_to_set_initstate=variable_names[variable_names %in% c("cc","KRAS","DSB")]
# # the values these nodes take on in different versions of model
initstate_values_list <- list(c(1,0,0), c(1,0,1), c(1,1,1))
## FIXED GENES (to what values these genes would be fixed: -1 is KO, +1 knock-in)
# genes are fixed by: setting initial state for 0/1, up/down rate for 0
fixed_genes=variable_names[variable_names %in% c("CHEK1","MAPKAPK2")]; indices_fixed_genes=list(c(0,0),c(-1,0),c(0,-1),c(-1,-1))

# simul from  "~/research/models/KRAS_DNA_repair_model/MaBoSS-Sampling/"
# input JSON file: 
sampling_settings_path = "~/research/models/KRAS_DNA_repair_model/MaBoSS-Sampling/examples/"
# results to: 
sampling_results_path="~/research/models/KRAS_DNA_repair_model/maboss_models/param_sampling_output"

# sampling_settings = fromJSON(paste(sampling_settings_path,"settings_example.json",sep=""))

col_df_all_conditions=c("name","value"); col_obj_states=c("name","value","type"); 
colnames_df_conditions_objectives=c("conditions","objectives")
# CONDITIONS (dataframe)
df_all_conditions=function_create_conditions_df_json(nodes_to_set_initstate, initstate_values_list,
                                                     fixed_genes,indices_fixed_genes)
sampling_settings_cell_lines=data.frame(matrix(NA,nrow=max(unique(df_all_conditions$cond_id)),ncol=1)); colnames(sampling_settings_cell_lines)="conditions"
sampling_settings_cell_lines$conditions=split(df_all_conditions[,colnames(df_all_conditions) %in% col_df_all_conditions], f=df_all_conditions$cond_id)

# PARAMETERS with values to sample (dataframe)
params_optim=c('$d_DSB', '$u_CDC25B'); colnames_optim_param_df=c("name", "values") # 
param_val_range = c(0.1,10^(-0.5),1,10^0.5,10,100); # this is assuming same paramvalues for all param.s
x=matrix(rep(param_val_range, length(params_optim)),ncol=length(params_optim), nrow=length(param_val_range))
sampling_settings_parameters=data.frame(matrix(NA,nrow=length(params_optim),ncol=length(colnames_optim_param_df))); colnames(sampling_settings_parameters)=colnames_optim_param_df
sampling_settings_parameters$name=params_optim; sampling_settings_parameters$values=split(x, rep(1:ncol(x), each = nrow(x)))

# JOIN and WRITE to JSON (list of dataframes)
list_sampling_settings_json=list(cell_lines = sampling_settings_cell_lines, parameters=sampling_settings_parameters)
write(toJSON(list_sampling_settings_json, pretty=T), paste(sampling_settings_path,"settings.json",sep=""))

# command for sampling
# MaBoSS-Sampling -c|--config CONF_FILE -s|--settings SETTINGS_FILE] [-o|--output RESULT_FILE] BOOLEAN_NETWORK_FILE

input_cfg="~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.cfg"
input_bnd="~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.bnd"
sampling_command = paste("./MaBoSS-Sampling -c ",input_cfg, 
                         " -s examples/settings.json -o ", sampling_results_path, "/sampling_results_json.csv ", input_bnd, sep = "")
# print command into SH file
write(sampling_command, "~/research/models/KRAS_DNA_repair_model/MaBoSS-Sampling/sampler_json.sh")

############################################################
# EXTRACT results

sampling_results=read_tsv(paste(sampling_results_path,"sampling_results_json.csv",sep = "/"))
nodes_data_constr = c("DSB", "CHEK1", "MAPKAPK2", "CDC25B", "cell_death")
if (max(sampling_results$`Condition Id`)<length(unique(sampling_results$`Condition Id`))) {
  sampling_results$`Condition Id`=sampling_results$`Condition Id`+1
}

id_cols=c("Condition Id","ParameterSet Id")
sampling_results_cut=function_sampling_results_cut(sampling_results,max_variables, nodes_data_constr,id_cols,"cond_id")
# sampling_results_cut_tidy=melt(sampling_results_cut,id.vars=c("cond_id","ParameterSet Id"))

# join exp vals to simul table
colnames_exp_values = c("variable","cond_id"); name_paramset_id="ParameterSet Id"; joining_cols=c("ParameterSet Id","rank_abs_dev","rank_sq_dev")
# DATAFRAME with errors
mean_errors_paramset = function_create_errors_df_parsampling(melt(sampling_results_cut,id.vars=c("cond_id","ParameterSet Id")),
                                exp_values, sampling_results_cut, colnames_exp_values,name_paramset_id)
# TIDY DATAFRAME with errors for plotting
n=2; errors=c("abs_dev", "sq_dev"); error_type_not_displayed=errors[c(1,2)[!c(1,2)==n]]; 
mean_errors_paramset_tidy = function_create_tidy_df_errors_parsampling(mean_errors_paramset, errors, error_type_not_displayed)
cols_display = c("rank_abs_dev", "rank_sq_dev")

####
# PLOT param values & MSE/MAE from data

# SAVE
height_width_vals=c(10,8); paramscan_plot_name="sampling_results_json"
postscript(paste(paramscan_path, paramscan_plot_name,"_",cols_display[n],".eps",sep=""), 
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
#
y_axis=cols_display[n]; cols_plot_df=c(y_axis,"variable","value","value_mod"); first_col_name="y_axis"
plot_df=function_create_plot_df_parsampling_errors(mean_errors_paramset_tidy,c(y_axis,"variable","value","value_mod"),"y_axis")
# demarc line
my_lines=function_my_lines("rank",y_axis,plot_df); theme1$axis.text.x=element_text(size=12); theme1$axis.text.y=element_blank()
ggplot(plot_df,aes(x=variable,y=y_axis,fill=value_mod)) + geom_tile(color="black",size=0.02) + 
  geom_text(aes(label=round(value,2)), size=4) + scale_fill_gradient2(low="green",mid="white",high="red") + 
  geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.5,inherit.aes=F,color="blue") + 
  ggtitle(paramscan_plot_name) + guides(fill=FALSE) + theme1
dev.off()

# PLOT errors w index of paramset
n=2; n_to_show=dim(mean_errors_paramset)[1]; mse=mean_errors_paramset[,errors[n]]; smallest_n=round(n_to_show*0.1)
title_mse=paste("mean",errors[n]); xlab_var="model index"
postscript(paste(paramscan_path, paramscan_plot_name,"_",cols_display[n],"_dots",".eps",sep=""), height=6,width=8,onefile=FALSE,horizontal=FALSE,paper='special')
function_mse_plot(mean_errors_paramset[,errors[n]], n_to_show,smallest_n,title_mse,xlab_var) # gsub("rank_","",y_axis)
dev.off()

# errors as function of pars with line plots
par_cols=colnames(mean_errors_paramset)[grepl("\\$",colnames(mean_errors_paramset))]
mean_errors_lineplot=melt(mean_errors_paramset[,colnames(mean_errors_paramset)[grepl("\\$",colnames(mean_errors_paramset)) | (colnames(mean_errors_paramset) %in% errors)]],
     measure.vars = errors)
r_cols=paste("r_",gsub(".*_","",par_cols),sep="")
for (i in 1:length(par_cols)) {
  mean_errors_lineplot[,r_cols[i]]=NA
  mean_errors_lineplot[,r_cols[i]]=mean_errors_lineplot[,par_cols[i]]/(1+mean_errors_lineplot[,par_cols[i]])
}

n=2
height_width_vals=c(6,12); font_size=16
postscript(paste(paramscan_path,"parscan_ERROR_gradients_",par_cols[n],".eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
char_str=par_cols[!par_cols %in% par_cols[n]]
ggplot(mean_errors_lineplot,
  aes(x=get(r_cols[n]),y=value,group=get(par_cols[!par_cols %in% par_cols[n]]),colour=factor(get(par_cols[!par_cols %in% par_cols[n]]))) ) +
  geom_line(size=1.5) + geom_point(size=2) + xlab(r_cols[n]) + facet_wrap(~variable,ncol=2,labeller=label_wrap_gen(multi_line=FALSE)) + ylab("error") +
   labs(col=char_str) + theme(strip.background = element_blank(),panel.grid=element_blank(),text=element_text(size=font_size),strip.text=element_text(size=font_size) ) 
dev.off()
# ,axis.text=element_text(size=font_size), legend.text=element_text(size=font_size),strip.text=element_text(size=font_size)

############################################################
# VISUALIZE value of NODES with LINE PLOTS as fcn of par. values

sampling_results_lineplot=melt(sampling_results[,colnames(sampling_results) %in% c(vars_to_plot,"Condition Id","ParameterSet Id") | 
                                                  grepl("\\$",colnames(sampling_results))], measure.vars = vars_to_plot)
# let's plot end value only (for CHEK1, MAPKAPK2), not max, because latter probably depends on par vals in complex ways 

par_cols=colnames(sampling_results_lineplot)[grepl("\\$",colnames(sampling_results_lineplot))]
r_cols=paste("r_",gsub(".*_","",par_cols),sep="")
for (i in 1:length(par_cols)) {
  sampling_results_lineplot[,r_cols[i]]=NA
  sampling_results_lineplot[,r_cols[i]]=sampling_results_lineplot[,par_cols[i]]/(1+sampling_results_lineplot[,par_cols[i]])
}

# as a function of $d_DSB
height_width_vals=c(10,12); max_time=200; font_size=16
n=1; char_str=par_cols[!par_cols %in% par_cols[n]]
postscript(paste(paramscan_path,"parscan_node_values_gradients_",par_cols[n],"_maxtime",max_time,".eps",sep=""),
             height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')

ggplot(sampling_results_lineplot, # [!sampling_results_lineplot$variable %in% c("CHEK1","DSB","MAPKAPK2"),]
  aes(x=get(r_cols[n]),y=value,group=interaction(variable,get(r_cols[!par_cols %in% par_cols[n]])),
      colour=variable, shape=factor(get(par_cols[!par_cols %in% par_cols[n]])))) +
  geom_line() + geom_point() + # scale_x_continuous(trans='log10') + 
  facet_wrap(~`Condition Id`,ncol=4,labeller=label_wrap_gen(multi_line=FALSE)) + xlab(r_cols[n]) + labs(shape=char_str) + 
  theme(strip.background = element_blank(),panel.grid=element_blank(),text=element_text(size=font_size),strip.text=element_text(size=font_size) )
dev.off()

  # +  legend(element_text())


############################################################

# VISUALIZE results of param sampling for 1 condition with HEATMAPS, all param sets
# read in
# filename_parscan="cellfate_30timestep_cond1.csv"; paramscan_path="~/research/models/KRAS_DNA_repair_model/maboss_models/param_sampling_output/"
# create table with mean values of each variable for each param set, with `ParameterSet Id` and the parameter value
cond_id_val=8
pscan_results_means=sampling_results_cut[sampling_results_cut$cond_id==cond_id_val, !colnames(sampling_results_cut) %in% "cond_id"]; pscan_results_means=pscan_results_means[,sort(colnames(pscan_results_means))]
# function_meanvalues_paramscan(paramscan_path,filename_parscan)
pscan_results_means_tidy=melt(pscan_results_means,measure.vars=colnames(pscan_results_means)[!colnames(pscan_results_means) %in% "ParameterSet Id"])

# SAVE PLOT
height_width_vals=c(14,10)
postscript(paste(paramscan_path,"parscan_node_values_cond",cond_id_val,".eps",sep=""),height=height_width_vals[1],width=height_width_vals[2],
                                       onefile=FALSE,horizontal=FALSE,paper='special')
# ggplot themes
###
###
theme=theme(axis.title=element_blank(),axis.text.x=element_text(size=14),axis.text.y=element_blank(),plot.title=element_text(hjust=0.5),
            legend.title=element_blank(),legend.text=element_text(size=8), axis.ticks.y=element_blank(), legend.box.background=element_rect(),
            legend.box.margin=margin(1,3,3,3), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
### different color code for parameters and variables
pscan_results_means_tidy$value_mod=pscan_results_means_tidy$value; pscan_results_means_tidy$value_mod[grepl("\\$",pscan_results_means_tidy$variable)]=rescale(log10(pscan_results_means_tidy$value[grepl("\\$",pscan_results_means_tidy$variable)]))*(-1)
display_vars=colnames(pscan_results_means)[colnames(pscan_results_means) %in% vars_to_plot | grepl("\\$",colnames(pscan_results_means)) | grepl("Parameter",colnames(pscan_results_means))]
# dividing line between params and vars
xval=length(unique(pscan_results_means_tidy$variable[grepl("\\$",pscan_results_means_tidy$variable)]))+0.5
my_lines=data.frame(x=xval,y=-0.5,xend=xval,yend=length(unique(pscan_results_means_tidy$`ParameterSet Id`))-0.5)
# HEATMAP PLOT: first columns parameters, columns on the right variables
ggplot(pscan_results_means_tidy[pscan_results_means_tidy$variable %in% display_vars,],aes(x=variable,y=`ParameterSet Id`,fill=value_mod)) +
  geom_tile(color="black",size=0.02) + scale_fill_gradient2(low="green",mid="white",high="red") +
  geom_text(aes(label=round(value,2)),size=4) + ggtitle(paste("nodes as function of paramvals, cond=",cond_id_val)) + 
  geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.5,inherit.aes=F,color="blue") + guides(fill=FALSE) + theme 
###
###
dev.off()

############################################################
# check to see if results are same as simulations above launched w separate MaBoSS runs
# simults from separate runs
all_simulation_maxpoints$cond_id = match(all_simulation_maxpoints$condition_folder, conditions_list)

# join two simul results
id_vars=c("cond_id", "variable", "value"); simul_colnames=c("sep_runs","joint_runs")
x=all_simulation_maxpoints[all_simulation_maxpoints$variable %in% vars_to_plot,id_vars]; x$simul_type=simul_colnames[1]
y=melt(sampling_results_cut[sampling_results_cut$`ParameterSet Id`==order(mean_errors_paramset$sq_dev)[4],], id.vars = id_vars[1], measure.vars = vars_to_plot)
y$simul_type=simul_colnames[2]; df_plot_simul=rbind(x,y)
# differences
df_plot_simul_differences=join(x,y, by=c("cond_id", "variable"))
colnames(df_plot_simul_differences)[colnames(df_plot_simul_differences) %in% "value"]=paste("value",simul_colnames,sep="_")
colnames(df_plot_simul_differences)[colnames(df_plot_simul_differences) %in% "simul_type"]=paste("simul_type",simul_colnames,sep="_")
df_plot_simul_differences$val_diffs=abs(df_plot_simul_differences[,paste("value",simul_colnames,sep="_")[1]]-df_plot_simul_differences[,paste("value",simul_colnames,sep="_")[2]])

# results next to each other
ggplot(df_plot_simul,aes(x=factor(cond_id),y=variable,fill=value)) + scale_fill_gradient(low="white",high="red",limits=c(0,1)) +
  geom_tile(color="black",size=0.1) + facet_wrap(~simul_type,ncol=length(unique(df_plot_simul$simul_type)),labeller=label_wrap_gen(multi_line=FALSE)) + geom_text(aes(label=round(value,1)), size=geom_text_size)
# DIFFERENCES
postscript(paste(paramscan_path,"joint_separate_runs_diffs.eps",sep=""),height=height_width_vals[1],width=height_width_vals[2],
                                       onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(df_plot_simul_differences, aes(x=factor(cond_id),y=variable, fill=val_diffs )) + scale_fill_gradient(low="white",high="red",limits=c(0,1)) +
  geom_tile(color="black",size=0.1) + geom_text(aes(label=round(val_diffs,1)), size=geom_text_size)
dev.off()

############################################################
############################################################
# OPTIMIZATION
# with simulated annealing code by Vincent, https://github.com/vincent-noel/MaBoSS-Optim
# file is in ~/research/models/KRAS_DNA_repair_model/MaBoSS-env-2.0/engine/sampling/MaBoSS-Optim
# command:
#  "./MaBoSS-Optim -c examples/Four_cycle.cfg -s *.json -o results.json *.bnd"

library(jsonlite)
# read in JSON file: reads into a list, with dataframes
json_df_results=fromJSON("~/research/models/KRAS_DNA_repair_model/MaBoSS-env-2.0/engine/sampling/MaBoSS-Optim/results.json")

# create JSON table with conditions for cell lines
# JSON table made up of 3 things: 
# conditions (fixed parameters and initial values) + objectives (values of variables we want to have) + optimization parameters

# DATAFRAME with CONDITIONS
df_all_conditions=function_create_conditions_df_json(nodes_to_set_initstate, initstate_values_list,fixed_genes,indices_fixed_genes)
col_df_all_conditions=c("name","value"); col_obj_states=c("name","value","type"); colnames_df_conditions_objectives=c("conditions","objectives")
# create DF with all conditions and objectives
max_variables=c("CHEK1","MAPKAPK2")
# DATAFRAME with CONDITIONS & OBJECTIVES
df_conditions_objectives=function_create_conditions_objectives_df(df_all_conditions,max_variables,exp_values,col_obj_states)

# OPTIMIZATION PARAMETERS (dataframe): colnames("name","initial","min","max","digits") # digits is precision of optimization
params_optim=c('$d_HRR2', '$d_DSB', '$u_CDC25B'); init_vals=1; paramval_range=c(1e-2,1e2);
colnames_optim_param_df = c("name", "initial", "min", "max", "digits"); digits_val=2
optim_param_df=function_create_optim_param_df(params_optim, init_vals, paramval_range, colnames_optim_param_df, digits_val)

# final JSON table: list with all conditions + objectives + params to optimize
list_json_optim = list(cell_lines=df_conditions_objectives, optimization_parameters=optim_param_df)

# write JSON with conditions + objectives + params to optimize
settings_json_path_name = "~/research/models/KRAS_DNA_repair_model/MaBoSS-Optim/examples/settings_test.json"
write(toJSON(list_json_optim, pretty=T), settings_json_path_name)
# fromJSON("param_scan/settings_test.json")

# launch optimisation
# ./MaBoSS-Optim -c *.cfg -s settings.json -o results.json *.bnd
input_cfg <- "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.cfg"
input_bnd <- "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.bnd"
json_results_path_name="~/research/models/KRAS_DNA_repair_model/MaBoSS-Optim/results/optim_results.json"
command_optim = paste("./MaBoSS-Optim -c", input_cfg, "-s", settings_json_path_name, "-o", json_results_path_name , input_bnd , sep=" ")

# ./MaBoSS-Optim -c ~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.cfg -s ~/research/models/KRAS_DNA_repair_model/maboss_models/param_scan/settings_test.json -o param_scan/optim_results.json ~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/krasmodel1/model1.bnd --lambda 10

############################################################
############################################################
# Linear fitting of MSE by parameter values: MSE = f(param_vals)
linear_fit = lm(sq_dev~`$d_DSB`+`$d_HRR2`+`$u_CDC25B`, data=mean_errors_paramset)
# linear_fit0 = lm(sq_dev~`$d_HRR2`, data=mean_errors_paramset); anova(linear_fit,linear_fit0)

# abs percentage error
# mean(abs(residuals(linear_fit)/mean_errors_paramset$sq_dev))
# mean abs error: mean(abs(residuals(linear_fit)))
plot(mean_errors_paramset$`ParameterSet Id`,abs(residuals(linear_fit)/mean_errors_paramset$sq_dev),type = "b")
# abs error
plot(mean_errors_paramset$`ParameterSet Id`,abs(residuals(linear_fit)),type = "b")
# it's periodic...?!
# 
# VALUES vs RESIDUALS
plot(mean_errors_paramset$`ParameterSet Id`,mean_errors_paramset$sq_dev, type="b", 
     ylim=c(min(residuals(linear_fit)),max(mean_errors_paramset$sq_dev)))
lines(mean_errors_paramset$`ParameterSet Id`,fitted(linear_fit),col="blue")
lines(mean_errors_paramset$`ParameterSet Id`,residuals(linear_fit),col="red")
lines(mean_errors_paramset$`ParameterSet Id`,rep(0,dim(mean_errors_paramset)[1]), type = 'l', lty=2)

# dataframe with means grouped by $d_HRR2
mean_errors_paramset %>% group_by(`$d_HRR2`) %>% summarize_all(mean)
# plot(log10(unique(mean_errors_paramset$`$u_CDC25B`)),
# as.vector(as.matrix((mean_errors_paramset %>% group_by(`$u_CDC25B`) %>% summarize_all(mean))[,4])))

############################################################
# Linear prediction of variable values from param vals
# correlations between outputs
truth_vector_vars=!grepl("\\$",colnames(pscan_results_means)) & !grepl("Parameter",colnames(pscan_results_means))
cor(pscan_results_means[,truth_vector_vars])
# # linear fitting
lm_fit = lm(as.matrix(pscan_results_means[,truth_vector_vars]) ~ ., pscan_results_means[,grepl("\\$",colnames(pscan_results_means))])
summary(lm_fit); plot(pscan_results_means$`ParameterSet Id`, lm_fit$residuals[,1]); Anova(lm_fit)
# # MSE
# 
# # variables are heavily correlated, do PCA
pscan_results_pca=prcomp(pscan_results_means[,truth_vector_vars],center=T,scale.= T); par_set_names=rownames(pscan_results_means[,truth_vector_vars])
# # plot(pscan_results_pca,type="l)"
ggbiplot(pscan_results_pca,ellipse=T) + # ,obs.scale=1,var.scale=1
  scale_color_discrete(name = '') + theme(legend.direction='horizontal',legend.position='top')

# example
# data(wine); wine.pca <- prcomp(wine,scale.=TRUE)
# ggbiplot(wine.pca, obs.scale=1, var.scale=1, groups=wine.class, ellipse=TRUE, circle=TRUE) +
#   scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top')
