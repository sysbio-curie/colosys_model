# started: 26/oct/2018
# Mihaly Koltai

# this is a script to integrate mutational/CNV data and cell line viability measurements into the model
# library("Rcpp", lib.loc="/bioinfo/users/mkoltai/R/x86_64-redhat-linux-gnu-library/3.4")
library(ggplot2)
library("stringr", lib.loc="/bioinfo/users/mkoltai/R/x86_64-redhat-linux-gnu-library/3.4")
library(reshape2); library(readr); library(tidyr); library(dplyr); 
library(pracma); library(gridExtra); library(nlme); library("pryr", lib.loc="/bioinfo/users/mkoltai/R/x86_64-redhat-linux-gnu-library/3.4")
setwd("~/research/models/KRAS_DNA_repair_model")

##############################
# 03/oct/2019
# we have received new data from Natalie/Charité, stored in ~/research/models/KRAS_DNA_repair_model/data
# the cell line data is originally in 'all Bliss Scores - new version - Kopie.xlsx', 
# I have manually pasted this into (roughly) tidy format, stored in "all_Bliss_Scores_tidydata.csv"

##############################
# extract data from all files of Dietlein paper

# viability data is in "/bioinfo/users/mkoltai/research/data/dietlein/mmc4/SupplFile2/DataSet2"
folderpath <- "~/research/data/cell_line_charite/" # "~/research/data/dietlein/mmc4/SupplFile2/DataSet2"
chk1i_conc <-c(0.00,0.1,0.25,0.5,1,2.50) #  c(5,1, 0.5, 0.2, 0.05, 0)
mk2i_conc <- c(0,0.5,1,2.50) # c(10,2, 0.5, 0.2, 0.05, 0)

################
# EXTRACT Dietlein data
# filelist<-list.files(folderpath,pattern = ".tsv") # ".txt"
# for (i in seq(1,length(filelist))) {
# viab_matr <- read_tsv(paste(folderpath,filelist[i],sep = "/"), col_names=F)
# colnames(viab_matr) <- chk1i_conc; viab_matr$mk2 <- mk2i_conc; viab_matr$cell_line <- str_replace(filelist[i],".tsv","")
# viab_matr_tidy <- melt(viab_matr,id.vars = c("cell_line","mk2"))
# colnames(viab_matr_tidy)[colnames(viab_matr_tidy) %in% "variable"] <- "chek1"
# if (i>1) {
#   viab_df_tidy_all <- rbind(viab_df_tidy_all, viab_matr_tidy)
# } else {
#   viab_df_tidy_all <- viab_matr_tidy
# }
# }
# viab_df_tidy_all$chek1 <- as.numeric(as.character(viab_df_tidy_all$chek1))
# write_csv(viab_df_tidy_all, paste(folderpath, "all_viability_data.csv",sep = "/"))
################

viab_df_tidy_all <- read_csv(paste(folderpath, "all_viability_data.csv",sep = "/") )

#####################################################################################
#####################################################################################
# read in new data from Natalie/Charite (Oct 2019)

# viab_df_orig=read_csv("data/viability/all_Bliss_Scores_tidydata.csv")
# 
# # we make column of MK2 concentrations numerical
# viab_df_orig$mk2i=as.numeric(gsub("uM","",viab_df_orig$mk2i))
# mk2i_conc=unique(viab_df_orig$mk2i)
# chk1i_conc <- as.numeric(unlist(strsplit(sapply(strsplit(colnames(viab_df_orig)[grepl("chk1i",colnames(viab_df_orig),ignore.case = T)],"_"), function(x) x[2]),"uM")))
# colnames(viab_df_orig)[grepl("chk1i",colnames(viab_df_orig),ignore.case = T)]=chk1i_conc
# # melt to tidy format
# viab_df_tidy_all <- melt(viab_df_orig, id=c("cell_line","mk2i","replicate","time")) 
# colnames(viab_df_tidy_all)[grepl("variable",colnames(viab_df_tidy_all))]="chk1i"
# viab_df_tidy_all=viab_df_tidy_all[,c("cell_line","replicate","time","mk2i","chk1i","value")]
# write_csv(viab_df_tidy_all,"data/viab_df_tidy_all.csv")

# READ IN if already exists
folderpath="data/"
viab_df_tidy_all_replicates=read_csv("data/viability/viab_df_tidy_all.csv")

#################

# take means of replicates
# colnames(viab_df_tidy_all)[!grepl("replicate|value",colnames(viab_df_tidy_all))]
viab_df_tidy_all_timepoints=viab_df_tidy_all_replicates %>% group_by(cell_line,time,mk2i,chk1i) %>% summarise(coeffvar=sd(value)/mean(value),value=mean(value))
viab_df_tidy_all_timepoints$cell_line=gsub("pool","\npool",viab_df_tidy_all_timepoints$cell_line)

# PLOT viability of ALL cell lines (all timepoints)
height_width_vals <- c(6,40); axis_text_x_size<-14; font_size=4.5 # col_num <- ncol(viab_df_tidy_all); 
# p="~/research/data/dietlein/"
postscript(paste(folderpath,"viability_data.eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
####
# all timepoints in one plot
ggplot(viab_df_tidy_all_timepoints, aes(x=factor(mk2i), y=factor(chk1i), fill=value)) + 
  facet_grid(time~cell_line) + # facet_wrap(cell_line, ncol=col_num, labeller=label_wrap_gen(multi_line=FALSE)) +
  scale_fill_gradient(low="red",high="white", limits=c(0,1), na.value="transparent") + 
  geom_tile(color="black",size=0.1) + geom_text(aes(label=round(value,1)), size=font_size) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text=element_text(size=axis_text_x_size),axis.title=element_text(size=20), strip.text=element_text(size=11)) + 
    xlab("MK2i") + ylab("CHEK1i")
dev.off()

# for one timepoint only
timepoint=24
height_width_vals <- c(8,16); axis_text_x_size<-14; font_size=4.5 # col_num <- ncol(viab_df_tidy_all_timepoints); 
postscript(paste(folderpath,"viability_data_t",timepoint,"h.eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(viab_df_tidy_all_timepoints[viab_df_tidy_all_timepoints$time==timepoint,], aes(x=factor(mk2i), y=factor(chk1i), fill=value)) + 
  facet_wrap(~cell_line, ncol=col_num, labeller=label_wrap_gen(multi_line=FALSE)) + # 
  scale_fill_gradient(low="red",high="white", limits=c(0,1), na.value="transparent") + 
  geom_tile(color="black",size=0.1) + geom_text(aes(label=round(value,2)), size=font_size) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       axis.text=element_text(size=axis_text_x_size),axis.title=element_text(size=20), strip.text=element_text(size=14)) + 
  xlab("MK2i") + ylab("CHEK1i")
dev.off()

##############################
# viab_df_tidy_all$cell_line=gsub("\n","",viab_df_tidy_all$cell_line)

# use one timepoint only
timepoint=72
viab_df_tidy_all=viab_df_tidy_all_timepoints[viab_df_tidy_all_timepoints$time==timepoint,!colnames(viab_df_tidy_all_timepoints) %in% c("coeffvar","time")]
colnames(viab_df_tidy_all)=c("cell_line","mk2","chek1","value")
# single agent curve (same for every cell line)
chek1_conc=as.vector(as.matrix(as.data.frame( unique(viab_df_tidy_all[viab_df_tidy_all$mk2==0,"chek1"]) ))) # rev( )
mk2_conc=as.vector(as.matrix(as.data.frame( unique(viab_df_tidy_all[viab_df_tidy_all$chek1==0,"mk2"]) ))) # rev( )
n_interpol <- 500
fine_chek1_conc <- seq(0,max(chek1_conc), max(chek1_conc)/n_interpol); fine_mk2_conc<-seq(0,max(mk2_conc),max(mk2_conc)/n_interpol)
mesh_fine <- meshgrid(fine_chek1_conc, fine_mk2_conc)

##############################
# CALCULATE ALL SYNERGIES

source("cell_line_data_model_integration_functions.R")
# many cell lines: output is dataframe for all cell lines w diff synergy values. one cell line:list with parameters
start_cell_line<-1; stop_cell_line<-length(unique(viab_df_tidy_all$cell_line)); rev_val=""
viab_df_tidy_all_synerg <- function_synergy_calculator(start_cell_line,stop_cell_line,viab_df_tidy_all,chek1_conc,mk2_conc, 
                            n_interpol, fine_chek1_conc, fine_mk2_conc, mesh_fine,rev_val)
# output
# list(viab_df_tidy_all_synerg_pred_mk2, viab_df_tidy_all_synerg_pred_chek1)

# area based synergies
# mk2 predicted for EC50 curve
viab_df_tidy_all_synergy_areas_pred_mk2 <- unique(viab_df_tidy_all_synerg[[1]][,c(3,
                                              which(grepl("synergy_area",colnames(viab_df_tidy_all_synerg[[1]]))))]) 
## chek1 predicted for EC50 curve
viab_df_tidy_all_synerg[[2]] <- viab_df_tidy_all_synerg[[2]][,colnames(viab_df_tidy_all_synerg[[1]])]
viab_df_tidy_all_synergy_areas_pred_chek1 <- unique(viab_df_tidy_all_synerg[[2]][,c(3,
                                              which(grepl("synergy_area",colnames(viab_df_tidy_all_synerg[[2]]))))])

# CREATE a TABLE with mean of local synergy measures and both pred_mk2 and pred_chek1 area-based synergy measures
all_synergy_measures_per_cell_line <- function_create_per_cell_line_table(viab_df_tidy_all_synerg)

folderpath=paste("data/synergy_",timepoint,"h/",sep="")
write_csv(all_synergy_measures_per_cell_line, paste(folderpath,"all_synergy_measures_per_cell_line.csv",sep=""))
write_csv(viab_df_tidy_all_synerg[[1]][,!grepl("synergy_area",colnames(viab_df_tidy_all_synerg[[1]]))],
          paste(folderpath,"all_synergy_measures_per_dose.csv",sep="") )

####################################################################
####################################################################
# SUMMARY PLOTS: about all calculations

# scatter plots comparing all synergy measures (local and global)
# axis limits
cutoff_val<-0.15
lims=t(as.data.frame(lapply(2:ncol(all_synergy_measures_per_cell_line), function(r) {
  c(quantile(all_synergy_measures_per_cell_line[,r],cutoff_val,na.rm=T),quantile(all_synergy_measures_per_cell_line[,r],1-cutoff_val,na.rm=T)) }) ))
rownames(lims)<-colnames(all_synergy_measures_per_cell_line)[2:ncol(all_synergy_measures_per_cell_line)]

list_comp_syn<-list( c(1,2), c(1,3), c(2,3), c(1,4), c(1,5), c(1,6), c(1,7), c(3,4), c(3,5), c(3,6), c(3,7) )
list_comp_syn_pred_chek1 <- list( c(1,2), c(1,3), c(2,3), c(1,8), c(1,9), c(1,10), c(1,11), c(3,8), c(3,9), c(3,10), c(3,11) )

# WHICH LIST!!!!!  
list_subplots = list_comp_syn_pred_chek1; # filename="~/research/data/dietlein/comparison_all_synergy_measures.eps"
filename <- function_plot_name(list_subplots, paste(folderpath,"comparison_all_synergy_measures.eps",sep=""))
mar_vals <- c(5,5,5,2); subplot_vals=c(3,4); cex_lab_val<-2.5;cex_main<-3;cex_axis_val=1.5;margin_params<-c(5,5,5,2)
postscript(filename,height=15,width=18,onefile=FALSE,horizontal=FALSE,paper='special')
function_comparison_all_synergy_measures(list_subplots,all_synergy_measures_per_cell_line,mar_vals,subplot_vals,cex_lab_val,cex_axis_val)
dev.off()

####################
# scatter plots of correlations between pred_CHK1 and pred_MK2 values of 4 area-based synergy measures
postscript(paste(folderpath,"scatterplot_synergy_predmk2_predchek1_comparison.eps",sep=""),
           height=15,width=18,onefile=FALSE,horizontal=FALSE,paper='special')
par(mfrow=c(2,2),mar=c(5,5,5,2)) # 
for (counter in 2:(ncol(viab_df_tidy_all_synergy_areas_pred_mk2)) ) {
  xlim_vals<-ylim_vals<-c(-5,10)
  if (counter==5) {xlim_vals<-ylim_vals<-c(-10,50)}
  title_str=paste(str_replace(colnames(viab_df_tidy_all_synergy_areas_pred_chek1)[counter],"synergy_area_",""), ", corr=",
  round(cor(viab_df_tidy_all_synergy_areas_pred_chek1[,counter],
            viab_df_tidy_all_synergy_areas_pred_mk2[,counter],use="pairwise.complete.obs"),2),sep = "")
  plot(viab_df_tidy_all_synergy_areas_pred_chek1[,counter], viab_df_tidy_all_synergy_areas_pred_mk2[,counter],
       xlab="CHK1 predicted", ylab="MK2 predicted", xlim=xlim_vals, ylim=ylim_vals, 
       cex=3, cex.lab=cex_lab_val,cex.axis=cex_axis_val); title(title_str,cex.main=4,font.main=1)
}
dev.off()

####################
# scatter plots of correlations between the 4 area-based synergy measures
list_cols<-list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4));
cex_val=3;cex_lab_val=2.5;cex_axis_val=1.5; cex_main_val=4;mtext_string="pred CHK1";mtext_val=2.5
#  for pred_CHK1
postscript(paste(folderpath,"scatterplot_synergy_area_measures_predchek1.eps",sep=""),
           height=15,width=18,onefile=FALSE,horizontal=FALSE,paper='special')
par(mfrow=c(3,2),mar=c(5,5,5,2))
function_synergy_area_compare_plot(viab_df_tidy_all_synergy_areas_pred_chek1,list_cols,cex_val,cex_lab_val,cex_axis_val,cex_main_val,mtext_string,mtext_val)
dev.off()
# for MK2
postscript(paste(folderpath,"scatterplot_synergy_area_measures_predmk2.eps",sep=""),
           height=15,width=18,onefile=FALSE,horizontal=FALSE,paper='special')
par(mfrow=c(3,2),mar=c(5,5,5,2)); mtext_string="pred MK2"
function_synergy_area_compare_plot(viab_df_tidy_all_synergy_areas_pred_mk2,list_cols,cex_val,cex_lab_val,cex_axis_val,cex_main_val,mtext_string,mtext_val)
dev.off()
####################

# HEATMAP of synergy values per cell line
# remove "outer" measure (because often v high values)
area_synergy_measures_per_cell_line=all_synergy_measures_per_cell_line[,c(1,
      which(!grepl("outer",colnames(all_synergy_measures_per_cell_line)) & grepl("synergy_area",colnames(all_synergy_measures_per_cell_line)) ) ) ]
colnames(area_synergy_measures_per_cell_line) <- str_replace(colnames(area_synergy_measures_per_cell_line),"synergy_area_","")
num_cols=sapply(area_synergy_measures_per_cell_line,is.numeric); rownames(area_synergy_measures_per_cell_line)<-c()

plotname=paste(folderpath,"area_synergy_measures.eps",sep="")
norm_methods=c("scale","max")

norm_method=""

# scale columns (mean-centering)
if (norm_method=="scale") {
  area_synergy_measures_per_cell_line[,num_cols] <- scale(area_synergy_measures_per_cell_line[,num_cols])
  plotname=paste(folderpath,"area_synergy_measures_scaled.eps",sep = "")
} else if (norm_method=="max") { # divide by max area
area_synergy_measures_per_cell_line[,num_cols]=area_synergy_measures_per_cell_line[,num_cols]/
  max(apply(abs(area_synergy_measures_per_cell_line[,num_cols]),2,max,na.rm = T))
plotname=paste(folderpath,"area_synergy_measures_maxdivided.eps",sep="")
} else {
  plotname=paste(folderpath,"area_synergy_measures.eps",sep="")
}

# HEATMAP with area-based synergy values
height_width_vals=c(16,8); axis_text_x_size=18
limit_vals=c(min(apply(area_synergy_measures_per_cell_line[,num_cols],2,min,na.rm = T)), 
             max(apply(area_synergy_measures_per_cell_line[,num_cols],2,max,na.rm = T)))
postscript(plotname, height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(melt(area_synergy_measures_per_cell_line[,c(1,which(grepl("chek1",colnames(area_synergy_measures_per_cell_line))))]), 
       aes(x=variable, y=factor(cell_line,levels=rev(unique(cell_line))), fill=value)) + 
  # scale_fill_gradient(low="white",high="red", limits=limit_vals, na.value="grey") + 
  scale_fill_gradient2(low="blue", mid="white", high="red",na.value="grey",midpoint = 0) + # , limits=c(-10,10), space = "Lab", 
  geom_tile(color="black",size=0.1) + geom_text(aes(label=round(value,1)), size=8) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x=element_text(size=axis_text_x_size,angle=90),axis.text.y=element_text(size=axis_text_x_size)) + xlab("") + ylab("")
dev.off()

####################################################################
####################################################################
# PLOT synergy area measures for ONE selected CELL LINE
#
# STRUCTURE of output
# list_pars_indiv_cell_line= list(l_pred_mk2, l_pred_chek1, ec50_chek1_mk2_values)
#
# l_pred_chek1 <- list(z,pred_ec50_mk2_restrictive, pred_ec50_chek1_restrictive, pred_ec50_mk2_finite, pred_ec50_chek1_finite_inner, 
#      pred_ec50_chek1_finite_outer, mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad)

library("purrr", lib.loc="/bioinfo/local/build/Centos/R/R-3.6.0/lib64/R/library")

# in a loop: n_list=c(4,7,9,11,17)
# for (n in n_list) { }

n=20; start_cell_line=n; stop_cell_line=n
# calculate curves and areas for one cell line
list_pars_indiv_cell_line <- function_synergy_calculator(start_cell_line,stop_cell_line,viab_df_tidy_all,chek1_conc,mk2_conc,
                                                         n_interpol,fine_chek1_conc,fine_mk2_conc,mesh_fine,rev_val)
# VALUES!!
l_pred_mk2=list_pars_indiv_cell_line[[1]]; l_pred_chek1=list_pars_indiv_cell_line[[2]]; ec50_chek1_mk2_values=list_pars_indiv_cell_line[[3]]

# PLOT 4 SYN MEASURES for ONE CELL LINE
file_path=paste(folderpath,"viability/synergy_",timepoint,"h/indiv_cell_lines_synergy_plots/",sep = "")
height_width_vals<-c(12,18); cex_lab_val<-2.5;cex_main<-3;cex_axis_val=1.5; margin_params<-c(5,5,5,2)
# library(pryr) # Rgraphviz
##########
# CHEK1 is predicted
postscript(paste(file_path,unique(viab_df_tidy_all$cell_line)[start_cell_line],"_pred_chek1_synergy_areas.eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# PLOTTING FUNCTION
xmax=max(l_pred_chek1[[length(l_pred_chek1)-1]])*1.1; ymax=3
function_synergy_area_plot_pred_chek1(height_width_vals,start_cell_line,xmax,ymax,
                                      l_pred_chek1,ec50_chek1_mk2_values,cex_lab_val,cex_main,cex_axis_val,margin_params)
dev.off()

##########
# MK2 predicted (usually worse)
# postscript(paste(file_path,unique(viab_df_tidy_all$cell_line)[start_cell_line],"_pred_mk2_synergy_areas.eps",sep=""),
#            height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# # PLOTTING FUNCTION
# xmax<-max(l_pred_mk2[[length(l_pred_mk2)-1]])*1.1; if (is.infinite(xmax)) {xmax=5.5}; ymax=3
# cex_lab_val=2.5; cex_main=3; cex_axis_val=1.5; margin_params=c(5,5,5,2)
# function_synergy_area_plot_pred_mk2(height_width_vals, start_cell_line,xmax,ymax,
#                                     l_pred_mk2, ec50_chek1_mk2_values,cex_lab_val,cex_main,cex_axis_val,margin_params)
# dev.off()


########################################################################################################################################
########################################################################################################################################
# COMPARE to synergy classification provided by Dietlein

cell_line_classification_paper <- read_tsv("~/research/data/dietlein/mmc2.tsv")
cell_line_classification_paper=cell_line_classification_paper[order(cell_line_classification_paper$`Cell Line Name`),]
# synergy 0 or 1
synergy_vector_paper=cell_line_classification_paper$`Chk1i/MK2isynergistic`[cell_line_classification_paper$`Cell Line Name` %in% 
                                                                              all_synergy_measures_per_cell_line$cell_line]
# scale by column or use raw data (better)
norm_method="scale"

if (norm_method=="scale") {
  df_synergies <- scale(all_synergy_measures_per_cell_line[,grepl("synergy_area",colnames(all_synergy_measures_per_cell_line))])
  plotname="~/research/data/dietlein/summary_plots/corr_paper_owncalcul_measures_scaled.eps"
  thresholds=c(-rev(10^(seq(-4,0.3,0.1))),0,10^(seq(-4,0.3,0.1))) # c(0,  )
} else {
  df_synergies <- all_synergy_measures_per_cell_line[,grepl("synergy_area",colnames(all_synergy_measures_per_cell_line))]
  plotname="~/research/data/dietlein/summary_plots/corr_paper_owncalcul_measures.eps"
  thresholds=10^(seq(-1,0.5,0.005))
}

# calculate correlations for different threshold values
l_corr_paper_calc = lapply(thresholds, function(x) { synergy_vector_calc <- (df_synergies > x); synergy_vector_calc[synergy_vector_calc]=1
cor(synergy_vector_calc,synergy_vector_paper,use="pairwise.complete.obs") } )
# put list into dataframe
corr_paper_calc=t(do.call(cbind.data.frame, l_corr_paper_calc)); rownames(corr_paper_calc)=thresholds; corr_paper_calc=corr_paper_calc[,c(1:5,7,6,8)]
# put into tidy DF, remove outer area measure
tidy_corr_paper_calc=melt(corr_paper_calc[,!grepl("outer",colnames(corr_paper_calc))])
# how many cell lines have a synergy (area) score? (not NA)
synergy_score_defined=nrow(df_synergies)-rowSums(apply(df_synergies, 1, is.na))
tidy_corr_paper_calc$value_non_na_cell_lines = synergy_score_defined[match(tidy_corr_paper_calc$Var2, names(synergy_score_defined))]
axis_text_x_size=12; title_str="raw values"; if (norm_method=="scale") { title_str="scaled" }
# PLOT
postscript(plotname, height=8,width=16,onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(tidy_corr_paper_calc, aes(x=Var1, y=value)) +
  facet_wrap(~Var2+value_non_na_cell_lines, ncol=length(unique(tidy_corr_paper_calc$Var2))/2, labeller=label_wrap_gen(multi_line=FALSE)) +
  geom_line() + xlab("threshold") + ylab("correlation") + # ylim(0.6,1) +
  geom_vline(xintercept=1, linetype = "dashed", colour = "red") +
  theme(text=element_text(size=20),axis.text=element_text(size=axis_text_x_size),strip.text.x=element_text(size=14),
        plot.title = element_text(hjust = 0.5)) + ggtitle(title_str) + if (norm_method!="scale") {scale_x_log10() }
dev.off()

# what are the cell lines that cause discrepancy?
if (norm_method=="scale") {
  best_thresholds_scaled <- apply(corr_paper_calc,2,which.max)
  binarized_synergystic_scaled=data.frame(matrix(NA, nrow=nrow(df_synergies),ncol=ncol(df_synergies)))
  colnames(binarized_synergystic_scaled)=colnames(df_synergies)
  for (i in 1:length(best_thresholds_scaled)) {
    truth_vector <- (df_synergies[,i] > thresholds[best_thresholds_scaled[i]]); truth_vector[truth_vector]=1
    binarized_synergystic_scaled[,i] <- truth_vector } # cor(synergy_vector_paper, binarized_synergystic_scaled,use="pairwise.complete.obs")
} else { 
  best_thresholds_raw <- apply(corr_paper_calc,2,which.max)
  binarized_synergystic_raw=data.frame(matrix(NA, nrow=nrow(df_synergies),ncol=ncol(df_synergies)))
  colnames(binarized_synergystic_raw)=colnames(df_synergies)
  for (i in 1:length(best_thresholds_raw)) {
    truth_vector <- (df_synergies[,i] > thresholds[best_thresholds_raw[i]]); truth_vector[truth_vector]=1
    binarized_synergystic_raw[,i] <- truth_vector } 
}

# HEATMAP of synergy values from Dietlein paper (1st column) and best-fit binarized versions of area-based measures from my calcul.s
height_width_vals=c(18,8); axis_text_x_size=10
plotname="~/research/data/dietlein/summary_plots/binary_heatmaps_best_corr_rawcalc.eps"
if (norm_method=="scale") { 
  plotname=str_replace(plotname,"rawcalc","scaled"); df_plot=binarized_synergystic_scaled } else {
    df_plot=binarized_synergystic_raw   }

df_plot$cell_line=all_synergy_measures_per_cell_line$cell_line
df_plot$syn_meas_paper=synergy_vector_paper; df_plot=df_plot[,!grepl("outer",colnames(df_plot))]
colnames(df_plot)=str_replace(str_replace(colnames(df_plot),"synergy_area_",""),"finite_","")
df_plot=df_plot[,c( ncol(df_plot)-1,ncol(df_plot),1:(ncol(df_plot)-2) )]
# PREPARE FILE
postscript(plotname, height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# CREATE PLOT
ggplot(melt(df_plot), aes(x=variable, y=factor(cell_line,levels=rev(unique(cell_line))), fill=value )) + 
  scale_fill_gradient2(low="white", high="red",na.value="grey") + # , limits=c(-10,10), space = "Lab", 
  geom_tile(color="black",size=0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                            axis.text.x=element_text(size=axis_text_x_size,angle=90)) + ylab("cell line") + 
  ggtitle(paste(title_str,"synergy score binarized"))
dev.off()

# HEATMAP with discrepancies only
df_plot_diff=df_plot[,-2];df_plot_diff[,2:ncol(df_plot_diff)]=df_plot_diff[,2:ncol(df_plot_diff)] - df_plot[,2]
plotname=str_replace(plotname,"binary_heatmaps_best","binary_heatmaps_discrep_best")
postscript(plotname, height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(melt(df_plot_diff), aes(x=variable, y=factor(cell_line,levels=rev(unique(cell_line))), fill=value )) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",na.value="grey",midpoint=0,limits=c(-1,1)) + # , limits=c(-10,10), space = "Lab", 
  geom_tile(color="black",size=0.1) + ylab("cell line") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(size=axis_text_x_size,angle=90), 
        plot.title = element_text(hjust = 0.5)) + ggtitle(paste(title_str,"synergy score difference (paper - own calculations)")) 
dev.off()

write_csv(df_plot, str_replace(plotname,".eps", ".csv"))

# which cell lines show discrepancy?
df_plot_diff[rowSums(abs(df_plot_diff[,2:ncol(df_plot_diff)])>0,na.rm = T)>0,]
probl_cell_lines=which(rowSums(abs(df_plot_diff[,2:ncol(df_plot_diff)])>0,na.rm = T)>0)

# plot these
file_path<-"~/research/data/dietlein/indiv_cell_lines_synergy_plots/problematic_cell_lines/"
if (norm_method=="scale") {file_path=str_replace(file_path,"problematic_cell_lines/","problematic_cell_lines/scaled/")}
function_single_cell_line_plotter(file_path, probl_cell_lines)

#####################
# PLOT
# xmax<-max(l_pred_mk2[[length(l_pred_mk2)-1]])*1.1; if (is.infinite(xmax)) { xmax=5.5}; ymax<-15
# plot(ec50_chek1_mk2_values$fine_chek1_conc,ec50_chek1_mk2_values$fine_mk2_conc,xlim=c(0,xmax),ylim=c(0,10.5))
# lines(l_pred_mk2[[2]],l_pred_mk2[[3]],col="red"); lines(l_pred_mk2[[4]],l_pred_mk2[[5]],col="blue");
# lines(l_pred_mk2[[4]],l_pred_mk2[[6]],col="green"); lines(l_pred_mk2[[7]],l_pred_mk2[[8]],col="orange")
#
# xmax <- max(l_pred_chek1[[length(l_pred_chek1)-1]])*1.1; ymax<-5.5;
# plot(ec50_chek1_mk2_values$fine_mk2_conc,ec50_chek1_mk2_values$fine_chek1_conc,xlim=c(0,xmax),ylim=c(0,ymax))
# lines(l_pred_chek1[[2]],l_pred_chek1[[3]],col="red"); lines(l_pred_chek1[[4]],l_pred_chek1[[5]],col="blue")
# lines(l_pred_chek1[[4]],l_pred_chek1[[6]],col="green"); lines(l_pred_chek1[[7]],l_pred_chek1[[8]],col="orange")

######################################################
# LOEWE: <1 synergy. Bliss: <0 synergy

# scatter plot of LOCAL synergy measures per DOSE (Bliss, Loewe)
library(gridExtra)
grid.arrange(
  ggplot(viab_df_tidy_all_synerg[[1]], aes(x=loewe_syn_linear, y=loewe_syn)) + geom_point() + geom_smooth(alpha=0.2) +
    xlim(0, 2) + ylim(0, 2),
  # bliss vs interpol loewe
  ggplot(viab_df_tidy_all_synerg[[1]], aes(x=bliss_syn, y=loewe_syn)) + geom_point() + xlim(-1/2,1/2) + ylim(0, 2) + geom_smooth(alpha=0.2),
  # size=2, shape=23
  # bliss vs linear loewe
  ggplot(viab_df_tidy_all_synerg[[1]], aes(x=bliss_syn,y=loewe_syn_linear))+geom_point()+xlim(-1/2,1/2)+ylim(0, 2)+geom_smooth(alpha=0.2),
  ncol=2)

