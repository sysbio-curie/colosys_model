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
sampling_settings_parameters=data.frame(matrix(NA,nrow=length(params_optim),ncol=length(colnames_optim_param_df))); 
colnames(sampling_settings_parameters)=colnames_optim_param_df
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
