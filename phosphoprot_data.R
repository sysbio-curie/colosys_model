library(ggplot2)
library("stringr", lib.loc="/bioinfo/users/mkoltai/R/x86_64-redhat-linux-gnu-library/3.4")
library(reshape2); library(readr); library(tidyr); library(dplyr); 
library(pracma); library(gridExtra); library(nlme); library("pryr", lib.loc="/bioinfo/users/mkoltai/R/x86_64-redhat-linux-gnu-library/3.4")

# script for older data in /data/users/mkoltai/research/data/phosphoproteomic/charite/plotting_script.R
# 
# # Phosphoprotein data from Charite, received 27/11/2018
# # (go down for newer data)
# # 
# phosprot_filepath = "~/research/data/phosphoproteomic/charite"
# sw620_wes_data=read_csv("~/research/data/phosphoproteomic/charite/run model SW620 1.11-shortened HDR 4.0 normalised to Vinculin.csv")
# colnames(sw620_wes_data)=str_replace(colnames(sw620_wes_data)," normalised to Vinculin","")
# # first number in sample column is CHEK1i, second is MK2i
# # columns with inhibitor concentrations
# sw620_wes_data$CHEK1_i <- as.vector(as.matrix(as.data.frame(strsplit(str_replace(sw620_wes_data$Sample,"SW620 ",""),"/"))[1,]))
# sw620_wes_data$MK2_i <- as.vector(as.matrix(as.data.frame(strsplit(str_replace(sw620_wes_data$Sample,"SW620 ",""),"/"))[2,]))
# sw620_wes_data=sw620_wes_data=sw620_wes_data[,!colnames(sw620_wes_data) %in% "Sample"]
# sw620_wes_data=sw620_wes_data[,c("CHEK1_i","MK2_i",colnames(sw620_wes_data)[!colnames(sw620_wes_data) %in% c("CHEK1_i","MK2_i")] )]
# sw620_wes_data[sw620_wes_data=="DMSO"]<-0
# 
# # take mean of first 2 rows or remove 2nd
# method="mean"
# if (method=="mean"){
#   sw620_wes_data_norm = sw620_wes_data %>% group_by(CHEK1_i,MK2_i) %>% summarise_all(funs(mean)) } else {
#     sw620_wes_data_norm = sw620_wes_data[!duplicated(sw620_wes_data[,grepl("_i",colnames(sw620_wes_data))]),]
#   }
# 
# # normalization by first or second row
# # x=as.vector(as.matrix(sw620_wes_data[1, !grepl("_i",colnames(sw620_wes_data))]))
# # sw620_wes_data[, !grepl("_i",colnames(sw620_wes_data))]=sweep(sw620_wes_data[, !grepl("_i",colnames(sw620_wes_data))],2,x,"/")
# # or divide by maximum value / subtract minimum + divide by range
# m=sw620_wes_data_norm[, !grepl("_i",colnames(sw620_wes_data_norm))]
# # sw620_wes_data_norm[, !grepl("_i",colnames(sw620_wes_data_norm))]=apply(m, 2, function(X) (X - min(X))/diff(range(X))) # function(X) X/max(X)
# 
# # sigmoidal normalization, as in Saez-Rodriguez Data-driven logical model of signal transduction (MSB 2009)
# m_min_div=sweep( abs(sweep(m,2, as.vector(as.matrix( apply(m,2,min) )), "-")), 2, as.vector(as.matrix(m[1,])), "/")
# n=2;
# # n=2; (m_min_div^n)/(m_min_div^n + 0.5^n)
# # k=5; 1/(1 + exp(-k*(m_min_div-0.5)))
# # max divided
# # (m/max(m) )/( m/max(m) + 0.05 )
# # B processed value
# B=( (m/max(m) )/( m/max(m) + 0.05 ) )*( (m_min_div^n)/(m_min_div^n + 0.5^n) )
# B[sweep(m,2, as.vector(as.matrix(m[1,])), "-")<0] = 1-B[sweep(m,2, as.vector(as.matrix(m[1,])), "-")<0]
# sw620_wes_data_norm[, !grepl("_i",colnames(sw620_wes_data_norm))]=B
# 
# # HEATMAP
# geom_text_size=5; colnum=4; height_width_vals=c(8,18)
# postscript(paste(phosprot_filepath, "/", "sw620_data_2018_nov_27", rep("_",sign(nchar(method))), method,".eps",sep=""),
#            height=height_width_vals[1],width=height_width_vals[2], onefile=FALSE, horizontal=FALSE, paper='special' )
# ggplot(melt(sw620_wes_data_norm, id.vars = colnames(sw620_wes_data_norm)[grepl("_i",colnames(sw620_wes_data_norm))] ),
#        aes(x=CHEK1_i, y=MK2_i, fill=value )) + scale_fill_gradient(low="white",high="red",limits=c(0,1)) +
#   geom_tile(color="black",size=0.1) + facet_wrap(~variable,ncol=colnum,labeller=label_wrap_gen(multi_line=FALSE)) +
#   geom_text(aes(label=round(value,1)), size=geom_text_size) # values on plot 
# # geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.25,inherit.aes=F,color="blue")
# dev.off()

############################################################
############################################################

# new data end of february 2019

# normalization was not done the proper way, as for many proteins vinculin values were normalized my <vinculin> of HT29 for other 2 cell lines too
# I corrected this in the excel file "all data - vinculin +norm_averages_only.xlsx" and exported it to
# "all data - vinculin +norm_averages_only_CORRECTED_NORMALISATION.csv"

phosprot_data_3cell_lines=read_csv(paste(phosprot_filepath,"all data - vinculin +norm_averages_only_CORRECTED_NORMALISATION.csv",sep = "/"))
# split condition column 
phosprot_data_3cell_lines=phosprot_data_3cell_lines %>% separate(condition,c("CHEK1i","MK2i"),"/")
# fill in empty
phosprot_data_3cell_lines[is.na(phosprot_data_3cell_lines)]=0; phosprot_data_3cell_lines[phosprot_data_3cell_lines=="DMSO"]=0

# normalization: each protein should be normalized by its highest value (Across the 3 cell lines?)
# remove plyr then reload dplyr!!!! (otherwise group_by doesnt work)
phosprot_data_3cell_lines = phosprot_data_3cell_lines[!phosprot_data_3cell_lines$CHEK1i %in% "etoposide",] %>% 
  group_by(variable) %>% mutate(maxdiv=Average/max(Average))

# index conditions: we need [0,0], [max,0], [0,max], [max,max]
phosprot_data_3cell_lines_rel_conds = phosprot_data_3cell_lines[(phosprot_data_3cell_lines$CHEK1i %in% "1") | 
    (phosprot_data_3cell_lines$MK2i %in% "2.5" & !(phosprot_data_3cell_lines$CHEK1i %in% "0.5")) |
    (phosprot_data_3cell_lines$CHEK1i %in% "0") & (phosprot_data_3cell_lines$MK2i %in% "0"),]
phosprot_data_3cell_lines_rel_conds$condition = as.integer(interaction(phosprot_data_3cell_lines_rel_conds$CHEK1i, phosprot_data_3cell_lines_rel_conds$MK2i))

# plot heatmap of results
postscript(paste(phosprot_filepath,"phosprot_data_heatmap_etopREMOVED_norm.eps",sep="/"),height=height_width_vals[1],
           width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
geom_text_size=18
theme1$axis.text.x=element_text(size=geom_text_size); theme1$axis.text.y=element_text(size=geom_text_size);  # my_lines
theme1$strip.text.x=element_text(size=geom_text_size*2)
ggplot(phosprot_data_3cell_lines_rel_conds, aes(x=factor(condition), y=factor(variable, levels=rev(unique(variable))), fill=maxdiv)) +
  scale_fill_gradient(low="white",high="red",limits=c(0,1)) +
  geom_tile(color="black",size=0.1) + facet_wrap(~cell_line,ncol=3,labeller=label_wrap_gen(multi_line=FALSE)) +
  geom_text(aes(label=round(maxdiv,2)), size=geom_text_size/2) + # values on plot
  theme1 # geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.25,inherit.aes=F,color="blue")
dev.off()

###################################################

# October/2019
# Aggregated phosphoprotein+FACS data from Charite
# stored in ~/data/wes_facs/summary_WES_data_and_cell_line_data.csv

setwd("~/research/models/KRAS_DNA_repair_model")
wes_facs_filepath="data/wes_facs/"
wes_facs_data=read_csv(paste(wes_facs_filepath,"summary_WES_data_and_cell_line_data.csv",sep = "") )

# clean column names
colnames(wes_facs_data)=gsub("\\(","",gsub(")","",gsub(" ","_",colnames(wes_facs_data))))
# remove 'µM' from treatment column
wes_facs_data$chek1i=wes_facs_data$treatment
wes_facs_data$mk2i=wes_facs_data$treatment
wes_facs_data$etoposide=wes_facs_data$treatment
# chek1 cleaning
wes_facs_data$chek1i[!grepl("Chk1i", wes_facs_data$chek1i)]=0
wes_facs_data$chek1i=as.numeric(sapply(strsplit(wes_facs_data$chek1i,"µM Chk1i"), '[[',1))
# mk2 cleaning
wes_facs_data$mk2i[!grepl("Mk2i", wes_facs_data$treatment)]=0
wes_facs_data$mk2i[grepl("Mk2i", wes_facs_data$mk2i)]=sapply(strsplit(wes_facs_data$mk2i[grepl("Mk2i", wes_facs_data$mk2i)],"µM Mk2i"), '[[',1)
wes_facs_data$mk2i[grepl("Chk1i", wes_facs_data$mk2i)]=sapply(strsplit(wes_facs_data$mk2i[grepl("Chk1i", wes_facs_data$mk2i)],"/ "),"[[",2)
wes_facs_data$mk2i=as.numeric(wes_facs_data$mk2i)
# etoposide cleaning
wes_facs_data$etoposide[!grepl("etoposide",wes_facs_data$etoposide)]=0
wes_facs_data$etoposide[grepl("etoposide",wes_facs_data$etoposide)]=sapply(strsplit(wes_facs_data$etoposide[grepl("etoposide",wes_facs_data$etoposide)],"µM"),"[[",1)
wes_facs_data$etoposide=as.numeric(wes_facs_data$etoposide)
# SAVE
write_csv(wes_facs_data,paste(wes_facs_filepath,"wes_facs_data.csv",sep=""))
# reorder columns
treatment_names = c("chek1i","mk2i","etoposide")
wes_facs_data=wes_facs_data[,c(which(colnames(wes_facs_data) %in% treatment_names),which(!colnames(wes_facs_data) %in% treatment_names))]
# MELT to tidy format
wes_facs_data_tidy <- melt(wes_facs_data, id=c(treatment_names,"treatment","cell_line")) 
# add a column indicating exper type
wes_facs_data_tidy$expertype=''; facs_truth_vals=grepl(pattern='[[:digit:]]h$',wes_facs_data_tidy$variable)
wes_facs_data_tidy$expertype[facs_truth_vals]='facs'; wes_facs_data_tidy$expertype[!facs_truth_vals]='wes'
# add a column indicating time (only for facs)
wes_facs_data_tidy$time[!facs_truth_vals]=''
wes_facs_data_tidy$time[facs_truth_vals] = str_match(as.character(wes_facs_data_tidy$variable[facs_truth_vals]),'\\d{2}h$')
wes_facs_data_tidy$time=as.numeric(gsub("h","",wes_facs_data_tidy$time))
# reorder columns
wes_facs_data_tidy=wes_facs_data_tidy[, c(which(!colnames(wes_facs_data_tidy) %in% "value"),which(colnames(wes_facs_data_tidy) %in% "value"))]
# save tidy data
write_csv(wes_facs_data_tidy,paste(wes_facs_filepath,"wes_facs_data_tidy.csv",sep=""))

# WES data: we need to normalize our values to [0,1] interval
# there are multiple ways to do it
# I try normalization from article of J Saez Rodr: https://cancerres.aacrjournals.org/content/suppl/2017/04/05/0008-5472.CAN-17-0078.DC1
# take the value for control for each P-protein --> for other conditions take log2 wrt control --> linearly rescale so as control is 0.5
# 
# Natalie already normalized the values by the DMSO (control) condition for each CELL LINE too, I'm not sure this is necessary,
# wrote email on 21/oct/2019 asking this. 
# For the moment I work with the data as it is
wes_facs_data_tidy = read_csv(paste(wes_facs_filepath,"wes_facs_data_tidy.csv",sep=""))

# take log2 wrt control, than rescale to [0,1]
wes_facs_data_tidy=wes_facs_data_tidy %>% group_by(cell_line,variable) %>% mutate(log2_scaled_value=(log2(value)-min(log2(value)))/(max(log2(value)) - min(log2(value)) ) )
# scaling without log2
wes_facs_data_tidy = wes_facs_data_tidy %>% group_by(cell_line,variable) %>% mutate(scaled_value=(value-min(value))/(max(value) - min(value) ) )
# scaling without etoposide values
wes_facs_data_tidy=wes_facs_data_tidy[wes_facs_data_tidy$etoposide==0,] %>% group_by(cell_line,variable) %>% mutate(scaled_value_no_etop=(value-min(value))/(max(value)-min(value) ) )
# scaling without etoposide values ACROSS cell lines
wes_facs_data_tidy=wes_facs_data_tidy[wes_facs_data_tidy$etoposide==0,] %>% group_by(variable) %>% mutate(scaled_value_no_etop_across_samples=(value-min(value))/(max(value)-min(value) ) )
# for FACS data this doesn't make sense, so delete it
scaled_value_names=colnames(wes_facs_data_tidy)[grepl("scaled",colnames(wes_facs_data_tidy))]
for (variable in scaled_value_names) { wes_facs_data_tidy[wes_facs_data_tidy$expertype %in% "facs",variable]=NA }

# truncated variable names
wes_facs_data_tidy$trunc_variable=sapply(strsplit(as.character(wes_facs_data_tidy$variable),"\\_"),'[[',1)
# enter 'etop' in mk2i/chk1i columns for etoposide condition
wes_facs_data_tidy$chek1i[wes_facs_data_tidy$etoposide>0]='etop'; wes_facs_data_tidy$mk2i[wes_facs_data_tidy$etoposide>0]='etop'

# save with scaled values as well
# wes_facs_data_tidy[,grepl("scaled",colnames(wes_facs_data_tidy))]=round(wes_facs_data_tidy[,grepl("scaled",colnames(wes_facs_data_tidy))],3)
# write_csv(wes_facs_data_tidy,paste(wes_facs_filepath,"wes_facs_data_tidy_scaled_vals.csv",sep=""))

wes_data_plot=wes_facs_data_tidy[wes_facs_data_tidy$expertype %in% "wes",]
wes_data_plot[,grepl("scaled",colnames(wes_data_plot))]=round(wes_data_plot[,grepl("scaled",colnames(wes_data_plot))],2)
wes_data_plot$diff_vals_norm = round(wes_data_plot$scaled_value - wes_data_plot$log2_scaled_value,2)

# PLOT WES DATA
height_width_vals=c(6,18); 
axis_text_x_size=18; geom_text_size=12 # geom_text_size=5; colnum=4; height_width_vals=c(8,18)
theme1=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title=element_text(size=18,hjust=0.5),
             axis.text=element_text(size=axis_text_x_size),axis.title=element_text(size=18), strip.text=element_text(size=18))
# theme1$axis.text.x=element_text(size=geom_text_size); theme1$axis.text.y=element_text(size=geom_text_size);  # my_lines
# theme1$strip.text.x=element_text(size=geom_text_size*2)
# with no etoposide or max etop
plot_var_name="scaled_value_no_etop_across_samples" # "scaled_value_no_etop" # scaled_value_no_etop scaled_value log2_scaled_value diff_vals_norm
min_max_vals=c(min(wes_data_plot[,plot_var_name]),max(wes_data_plot[,plot_var_name]))
if (min_max_vals[1]>=0) {min_max_col=c("white","red")} else {min_max_col=c("blue","white")}
postscript(paste(wes_facs_filepath,"phosprot_data_heatmap_",plot_var_name,".eps",sep=""),height=height_width_vals[1],
           width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(wes_data_plot, aes_string(x="mk2i",y="chek1i",fill=plot_var_name)) +
  scale_fill_gradient(low=min_max_col[1],high=min_max_col[2],limits=min_max_vals,guide=FALSE) + 
  geom_tile(color="black",size=0.1) + facet_grid(cell_line~trunc_variable) + 
  geom_text(aes_string(label=plot_var_name),size=geom_text_size/2) + xlab("MK2i") + ylab("CHEK1i") + # values on plot
  theme1
dev.off()

################################
# PLOT FACS data
facs_data_plot=wes_facs_data_tidy[wes_facs_data_tidy$expertype %in% "facs",]
facs_data_plot$trunc_variable=unlist(strsplit(as.character(facs_data_plot$variable),'_\\d{2}h$'))
height_width_vals=c(6,8)

# save separately for each variable
plot_list=list()
for (k in 1:length(unique(facs_data_plot$trunc_variable))) {
val_name=unique(facs_data_plot$trunc_variable)[k]
# postscript(paste(wes_facs_filepath,"FACS_data_heatmap_",val_name,".eps",sep=""),height=height_width_vals[1],
#            width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
min_max_vals=c(0,max(facs_data_plot[facs_data_plot$trunc_variable %in% val_name,"value"],na.rm = T)); 
p=ggplot(facs_data_plot[facs_data_plot$trunc_variable %in% val_name,], aes_string(x="mk2i",y="chek1i",fill="value")) +
  scale_fill_gradient(low="white",high="red",limits=min_max_vals,guide=FALSE) + 
  geom_tile(color="black",size=0.1) + facet_grid(cell_line~time) + ggtitle(val_name) + 
  geom_text(aes_string(label="value"),size=geom_text_size/2) + xlab("MK2i") + ylab("CHEK1i") + # values on plot
  theme1
plot_list[[k]] = p
# dev.off()
}

for (k in 1:length(unique(facs_data_plot$trunc_variable))) {
  val_name=unique(facs_data_plot$trunc_variable)[k]
  postscript(paste(wes_facs_filepath,"FACS_data_heatmap_",val_name,".eps",sep=""),height=height_width_vals[1],
             width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
  print(plot_list[[k]])
  dev.off()
}

################################
# Mutation data

mutation_data=read_csv("data/panel_seq_data.nbu_27092019.csv")
colnames(mutation_data)=gsub(" ","_",colnames(mutation_data))
colnames(mutation_data)[grepl("MSI",colnames(mutation_data))]="MSI_MSS"
mutation_data=mutation_data[,!colnames(mutation_data) %in% "Media"]

# integrate other mutations as well
other_mutations=lapply(strsplit(mutation_data$other_mutations_CRC_panel,"; "), function(x) {strsplit(x," ")}) 
other_mutations_table=data.frame(do.call(rbind, 
  lapply(1:length(other_mutations), function(x) {cbind(rep(x,length(other_mutations[[x]])),t(matrix(unlist(other_mutations[[x]]),nrow=2)))})),stringsAsFactors=F)
colnames(other_mutations_table)=c("Cell_line","variable","value")
other_mutations_table=other_mutations_table[!is.na(other_mutations_table$variable),]
other_mutations_table$Cell_line=mutation_data$Cell_line[as.numeric(other_mutations_table$Cell_line)]
other_mutations_table$other_mutations_CRC_panel = mutation_data$other_mutations_CRC_panel[match(other_mutations_table$Cell_line,mutation_data$Cell_line)]

# melt original table into tidy format
mutation_data_tidy=melt(mutation_data, id=c("Cell_line","other_mutations_CRC_panel")) 

# merge two tables
if (sum(mutation_data_tidy$variable %in% "APC")==0) {
  mutation_data_tidy=rbind(mutation_data_tidy[,c("Cell_line","variable","value","other_mutations_CRC_panel")], other_mutations_table)  
  # order by cell lines
  mutation_data_tidy=mutation_data_tidy[order(mutation_data_tidy$Cell_line,mutation_data_tidy$variable),]
  mutation_data_tidy=mutation_data_tidy %>% group_by(Cell_line,variable) %>% summarize(value=paste(value,collapse=",\n "),num_mut=n())
  mutation_data_tidy$num_mut[mutation_data_tidy$value %in% c("wt","?","MSS","NA")]=""
}

mutation_data_tidy$color=NA
mutation_data_tidy$color[is.na(mutation_data_tidy$value)| mutation_data_tidy$value %in% "NA" | mutation_data_tidy$value %in% "?"]="no data"
mutation_data_tidy$color[mutation_data_tidy$value %in% c("wt","MSS")]="wt or MSS"
mutation_data_tidy$color[is.na(mutation_data_tidy$color)]=paste("mutation", mutation_data_tidy$num_mut[is.na(mutation_data_tidy$color)])
mutation_data_tidy$color=factor(mutation_data_tidy$color)

# SAVE
write_csv(mutation_data_tidy,"data/mutation/mutation_data_tidy.csv")

# PLOT
axis_text_x_size=14; geom_text_size=12 # geom_text_size=5; colnum=4; 
height_width_vals=c(8,18)
theme1=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title=element_text(size=22,hjust=0.5),
        axis.text=element_text(size=axis_text_x_size), axis.text.x=element_text(angle=90), 
        axis.title=element_text(size=18), strip.text=element_text(size=18),
        legend.text=element_text(size=18))
#
postscript("data/mutation/mutation_data.eps",height=height_width_vals[1],
           width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(mutation_data_tidy, aes(x=variable,y=Cell_line,fill=color)) +
  scale_fill_manual(values=c("red1","red2","red3","red4", "grey", "green"),name="") + 
  geom_tile(color="black",size=0.1) + ggtitle("mutational data") + 
  geom_text(aes(label=num_mut),size=6) + xlab("genes or MSI") + ylab("cell lines") + 
  theme1
dev.off()

### 
# most common mutations
mutation_num_genes=melt(mutation_data_tidy %>% group_by(variable) %>% 
                          summarize(mult_count=sum(as.numeric(num_mut),na.rm=T), single_count=sum(as.numeric(num_mut)>0,na.rm=T)))
colnames(mutation_num_genes)[2]="count_type"

n_cell_lines=length(unique(mutation_data_tidy$Cell_line))
filter_truthvals=!is.na(mutation_num_genes$value) & 
  mutation_num_genes$variable %in% unique(mutation_num_genes$variable[mutation_num_genes$value>2])

postscript("data/mutation/mutation_data_occurrence.eps",height=8,width=8,onefile=FALSE,horizontal=FALSE,paper='special')
ggplot(mutation_num_genes[filter_truthvals,],aes(x=variable,y=value/n_cell_lines,fill=count_type)) + 
  geom_bar(stat="identity", color="black", position="dodge") + coord_flip() + 
  ylab("% occurrence") + xlab("") + scale_y_continuous(breaks=seq(0,1,0.1)) + theme1 + 
  theme(legend.position = c(0.8, 0.8), legend.text = element_text("sdsd")) +
  geom_hline(yintercept=0.2,linetype="dashed",color="blue",size=1.6)
dev.off()

write_csv(mutation_num_genes, "data/mutation/mutation_num_genes.csv")

unique(mutation_num_genes$variable[mutation_num_genes$value>4])