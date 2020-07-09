# functions for the script in maboss_wrapper_plotting.R

#######################################################################################
#######################################################################################
# FUNCTIONS

# function to create a folder for a model (.bnd file) with different cfg files (KOs/K-ins and initial states)
function_cfg_file_creator <- function(input_bnd_file, input_cfg_file, filename_tag_list, nodes_to_set_initstate, 
                                      initstate_values_list, fixed_genes, indices_fixed_genes,copy_move){

  file_counter <- 0; cfg_list <- c()
  # generate different filenames INIT STATE variants for diff. mutants (DNA damage)
  for (is_counter in 1:length(filename_tag_list)) {
    output_nametag <- filename_tag_list[is_counter]; initstate_values <- initstate_values_list[[is_counter]]
    # SET INIT STATES
    input_cfg_file_mutants_perturbs <- function_add_initial_states_rewrite(input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag)
    file_counter <- file_counter+1; cfg_list[file_counter] <- input_cfg_file_mutants_perturbs
    
    if (length(indices_fixed_genes)>0) {
    for (counter in 1:length(indices_fixed_genes) ) {
      OFF_mutants <- fixed_genes[indices_fixed_genes[[counter]]<0]
      ON_mutants <- fixed_genes[indices_fixed_genes[[counter]]>0]
      if ( length(ON_mutants)==0 ) {filetag <- paste("_fixed",paste0(OFF_mutants, "off",collapse = "_" ), sep="_")} 
      if ( length(OFF_mutants)==0 ) {filetag <- paste("_fixed",paste0(ON_mutants, "on",collapse = "_" ), sep="_")} 
      if ( length(OFF_mutants)>0 && length(ON_mutants)>0 ) {
        filetag <- paste("_fixed", paste0( ON_mutants, "on",collapse = "_"), paste0(OFF_mutants, "off",collapse = "_" ), sep="_")
      }
      # generate mutant CFG
      filename <- function_add_mutant_nodes(input_cfg_file_mutants_perturbs, ON_mutants, OFF_mutants, filetag)
      # store filenames
      file_counter <- file_counter+1; cfg_list[file_counter] <- filename
    }
    }
  }

  # move into directory, create cfg files, this is optional, only if argument "copy_move" non-empty
  if (nchar(copy_move)>0){
  dir_name <- gsub(".cfg","",input_cfg_file); 
  if (dir.exists(dir_name)==F) {
    dir.create(dir_name) 
  } else {
    dir_name <- paste(dir_name, "_1", sep = ""); dir.create(dir_name) 
  }
  file.rename(cfg_list, paste( dir_name, cfg_list, sep = "/") )
  system(paste("cp", input_bnd_file, dir_name))
  }
  
}

#######################################################################################

simulator_function <- function(input_cfg_file, input_bnd_file, dir_name, go_maboss_path){
  
  # dir_name <- gsub(".cfg","",input_cfg_file); 
  sh_command_vars <- paste(go_maboss_path, "\n source MaBoSS.env \n cd ", paste( getwd(), dir_name, sep = "/"), sep="")
  # generate list of commands
  # cfg_list <- list.files(dir_name, pattern = ".cfg")
  sh_command <- paste(sh_command_vars, paste0(paste("time MBSS_FormatTable.pl", input_bnd_file, input_cfg_file ), collapse = "\n"), "cd ..", sep = "\n")
  # write shell script with list of commands in folder
  write(sh_command, paste(dir_name, "run_simul.sh",sep="/"))
  # RUN simulations
  system(paste("cd", dir_name, "\n", "sh run_simul.sh") )
  
}

#######################################################################################

function_add_mutant_nodes <- function(input_cfg_fileName, ON_mutants, OFF_mutants, filename_tag) {
  
  file_as_character_string <- readChar(input_cfg_fileName, file.info(input_cfg_fileName)$size)
  output_file_name <- paste( gsub(".cfg","",input_cfg_fileName),  filename_tag, ".cfg", sep="")
  
  if (ON_mutants!="" && length(ON_mutants)>0){
    for (counter in 1:length(ON_mutants)){
      file_as_character_string <- gsub( paste("\\$d_", ON_mutants[counter], "=1;", sep=""),
                                        paste("\\$d_", ON_mutants[counter], "=0;", sep=""), file_as_character_string )
      if (grepl( paste(ON_mutants[counter], ".istate=.*;", sep = ""), file_as_character_string)){
        file_as_character_string <- gsub( paste(ON_mutants[counter] ,"[.]istate=\\d;", sep=""),
                                          paste(ON_mutants[counter] ,".istate=1;", sep=""), file_as_character_string)
      }
    }
  }
  
  if (OFF_mutants!="" && length(OFF_mutants)>0){
    for (counter in 1:length(OFF_mutants)){
      file_as_character_string <- gsub( paste("\\$u_", OFF_mutants[counter], "=1;", sep=""),
                                        paste("\\$u_", OFF_mutants[counter], "=0;", sep=""), file_as_character_string )
      if (grepl( paste(OFF_mutants[counter], ".istate=.*;", sep = ""), file_as_character_string)){
        file_as_character_string <- gsub( paste(OFF_mutants[counter] ,"[.]istate=\\d;", sep=""),
                                          paste(OFF_mutants[counter] ,".istate=0;", sep=""), file_as_character_string)
      }
    }
  }
  fileConn <- file(output_file_name); write(file_as_character_string, fileConn); close(fileConn)
  output_file_name
}

#######################################################################################

function_change_rates <- function(input_cfg_files, node_list, up_rates, down_rates) {
  
  for (input_cfg_fileName in input_cfg_files) {

  file_as_character_string <- readChar(input_cfg_fileName, file.info(input_cfg_fileName)$size)
  output_file_name <- input_cfg_fileName # paste( gsub(".cfg","",input_cfg_fileName),  filename_tag, ".cfg", sep="")

    for (counter in 1:length(node_list)){
      
      if (!is.na(down_rates[counter])) {
        file_as_character_string <- gsub(paste("\\$d_", node_list[counter], "=\\d+;", sep=""), 
                                      paste("\\$d_", node_list[counter], "=",down_rates[counter],";", sep=""), file_as_character_string)
      }
      if (!is.na(up_rates[counter])) {
        file_as_character_string <- gsub(paste("\\$u_", node_list[counter], "=\\d+;", sep=""), 
                                      paste("\\$u_", node_list[counter], "=",up_rates[counter],";", sep=""), file_as_character_string)
      }
    }

  fileConn <- file(output_file_name); write(file_as_character_string, fileConn); close(fileConn)
  input_cfg_fileName
  }
}

# gsub(paste("\\$d_", node_list[counter], "=\\d+;", sep=""), 
# paste("\\$d_", node_list[counter], "=",down_rates[counter],";", sep=""), "$d_dna_dam=5;")

#######################################################################################

function_add_initial_states_rewrite <- function(input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag){
  # read in cfg file that will be rewritten
  input_cfg_file_char <- unlist(strsplit(readChar(input_cfg_file, file.info(input_cfg_file)$size), "\n"))
  # rows with init states
  init_state_nodes <- input_cfg_file_char[grepl("istate", input_cfg_file_char )]
  # gsub("=.*","", nodes)
  nodes_to_replace <- init_state_nodes[gsub("\\..*","", init_state_nodes) %in% nodes_to_set_initstate]
  # if order of nodes different in .cfg file, need to reorder initstate values accordingly
  if (identical(order(nodes_to_replace), order(nodes_to_set_initstate))==F){
    initstate_values <- initstate_values[match(gsub("\\..*","", nodes_to_replace), nodes_to_set_initstate)]
  }
  
  init_state_nodes[gsub("\\..*", "", init_state_nodes) %in% nodes_to_set_initstate] <- unlist(lapply(1:length(initstate_values), 
                                           function(r) {gsub("=.*;", paste("=",initstate_values[r],";",sep = ""), nodes_to_replace[r]) } ))
  # replace rows of to be modified initstates with new values
  input_cfg_file_char[grepl("istate", input_cfg_file_char )] <- init_state_nodes
  
  output_file_name_initstates <- paste( gsub(".cfg","",input_cfg_file), output_nametag, ".cfg", sep = "")
  system(paste("cp ", input_cfg_file, output_file_name_initstates))
  # paste("\n\n", paste0(init_nodes, ".istate=", init_values, ";",collapse="\n"), sep = "")
  write(input_cfg_file_char , output_file_name_initstates )
  output_file_name_initstates
}

#######################################################################################

# with this function we can set the init. states of selected nodes AND that of all other nodes to OFF/ON/random
# (apart from those we specifically set) to all OFF, ON or random

function_add_initial_states <- function(input_cfg_file, nodes_to_set_initstate, initstate_values, all_OFF_ON_random, output_nametag){
  if (any(grepl(all_OFF_ON_random, c("OFF","ON", "on", "off","random")))) {

    # get list of nodes by finding first and last row with "$" sign
    input_cfg_file_char <- readChar(input_cfg_file, file.info(input_cfg_file)$size)
    nodes <- unlist(strsplit(input_cfg_file_char,"\n"))[grepl("\\$",unlist(strsplit(input_cfg_file_char,"\n")) )]
    nodes <- unique(gsub("\\$u_","",gsub("\\$d_","",gsub("=.*","", nodes))))

    if (any(grepl(all_OFF_ON_random, c("ON", "on")))  ){
      init_values <-rep(1, length(nodes))
      init_values[nodes %in% nodes_to_set_initstate] <- initstate_values
      init_nodes <- as.character(nodes)
    }
    if (any(grepl(all_OFF_ON_random, c("OFF", "off")) ) ){
      init_values <-rep(0, length(nodes))
      init_values[nodes %in% nodes_to_set_initstate] <- initstate_values
      init_nodes <- as.character(nodes)
    }
    if (all_OFF_ON_random=="random"){
      init_values <- initstate_values
      init_nodes <- nodes_to_set_initstate
    }
    
    output_file_name_initstates <- paste( gsub(".cfg","",input_cfg_file), output_nametag, ".cfg", sep = "")
    
    if (nchar(output_nametag)>0){
    system(paste("cp ", input_cfg_file, output_file_name_initstates))
    }
    # avoid duplicates, are init states already specified?
    init_nodes_test <- paste0(init_nodes, ".istate=", init_values)
    truth_vals <- sapply(lapply(1:length(init_nodes_test), 
                          function(x) {grepl( init_nodes_test[x], unlist(strsplit(input_cfg_file_char,"\n")) )} ), sum)>0
    
    initstates_string <- paste0(paste0(init_nodes, ".istate=", init_values, ";")[!truth_vals],collapse="\n")
    # initstates_string <- paste0(init_nodes, ".istate=", init_values, ";",collapse="\n")

    write( paste("\n", initstates_string, sep = ""), output_file_name_initstates, append = T)
    output_file_name_initstates
  }
  else {
    print("ERROR! specify init states for other nodes: on / off / random")
  }
}
#######################################################################################

function_internal_visible_nodes <- function(input_cfg_file, visible_nodes){
  
  input_cfg_file_char <- readChar(input_cfg_file, file.info(input_cfg_file)$size)
  split_file <- unlist(strsplit(input_cfg_file_char,"\n"))
  internal_nodes <- split_file[grepl("is_internal",split_file )]
  nodes_to_modif <- sapply(lapply(visible_nodes, function(n) grepl(n,internal_nodes) ), which)
  if ( length(nodes_to_modif)>0 ){
    split_file[grepl("is_internal",split_file )][nodes_to_modif] <- gsub("TRUE","FALSE",internal_nodes[nodes_to_modif])
    write(split_file, input_cfg_file)
  }
}

#######################################################################################

# saving phenotype plots

plot_dynamics_df <- function(df,foldername, plot_levels, legend_font_size, plot_width,plot_height,line_thickness) {

  plot_levels <- factor(sort(colnames(df)[2:ncol(df)]))
  if (any(grepl("nil",as.character(plot_levels)))) {
    plot_levels <- plot_levels[c(plot_levels[-which(plot_levels=="nil")], plot_levels[which(plot_levels=="nil")])]
  }
  
  tidy_df <- melt( df, id=c("Time"))
  
  if (length(plot_levels)<=9){
    require(RColorBrewer)
    myColors <- brewer.pal(length(plot_levels), "Set1")
    names(myColors) <- plot_levels
    colScale <- scale_colour_manual(name = "variable",values = myColors)
  }
  temporal <- ggplot( tidy_df, aes(x=Time,y=value, color=variable)) + 
    geom_line(size=line_thickness) + ylab("Phenotype probablity") +
    theme( legend.position = "right",legend.title = element_blank(), text=element_text(size=legend_font_size) ) + 
    guides(colour = guide_legend(nrow = length(unique(tidy_df$variable)))) + 
    scale_colour_discrete(drop=TRUE, limits = levels(plot_levels)) 
  if (length(plot_levels)<=9){
    temporal <- temporal+ colScale
  }
  ggsave(filename = paste(foldername, paste(deparse(substitute(df)),"_time_evolution.png", sep = ""), sep = "/"), plot=temporal, width = plot_width, height = plot_height)
  # return(temporal)
  
  barplot_data <- as.data.frame(t( df[ nrow(df), which(df[nrow(df), ] > 0)[which(df[nrow(df), ] > 0)>1]]))
  colnames(barplot_data)[1] <- "value"
  barplot_data$ord <- rownames(barplot_data)
  levels(barplot_data$ord ) <- levels(plot_levels)
  # barplot_data
  barplot_elementary_states <- ggplot(barplot_data, aes(x = ord, y = value, label = value*100, fill=ord ) ) +
    geom_col() +  geom_text(size = legend_font_size/3, position = position_stack(vjust = 0.5)) + coord_flip() + 
    theme(axis.ticks=element_blank(), panel.grid=element_blank(),  # , axis.title= element_blank()
          legend.position = "none",legend.title = element_blank(), text=element_text(size=legend_font_size) ) + ylab("Phenotype probablity") + xlab("") + 
    scale_fill_discrete(drop=TRUE, limits = levels(plot_levels) )
  
  ggsave(filename = paste(foldername, paste(deparse(substitute(df)),"_barplot.png", sep = ""), sep = "/"), 
         plot=barplot_elementary_states, width=plot_width, height=plot_height) 
}

#######################################################################################

multiple_folder_plotter <- function(dir_name, output_folders, plot_width, plot_height, plotting_threshold,line_thickness,font_size ){
  
  setwd(dir_name)
  
  for (counter in 1:length(output_folders) ) {
    output_folder <- output_folders[counter]
    output_filename <- paste(output_folder, gsub("^\\.","", output_folder), "_probtraj_table.csv", sep = "")
    
    probtraj_table <- read.csv(output_filename, sep = "\t", stringsAsFactors = F,check.names=FALSE)
    if ( all(is.na(probtraj_table[,ncol(probtraj_table) ])) ) {
      probtraj_table <- probtraj_table[,-ncol(probtraj_table) ]
    }
    probtraj_table_only_probs <- probtraj_table[, c(1:4,  seq(5, ncol(probtraj_table), by=2) )]
    colnames(probtraj_table_only_probs) <- gsub("\\]","",gsub("Prob\\[","", colnames(probtraj_table_only_probs) ) )
    
    probtraj_table_only_probs <- probtraj_table_only_probs[,c(1,5:ncol(probtraj_table_only_probs))]

    # display all phenotype trajectories
    probtraj_table_only_probs_TH <- create_plots_probtraj_table(output_folder, plot_width, plot_height, plotting_threshold,line_thickness,font_size)
    
    elementary_states_dynamics <- from_phenotypes_to_state_probabilities( probtraj_table_only_probs )
    
    plot_levels <- factor(sort(colnames(elementary_states_dynamics)[2:ncol(elementary_states_dynamics)]))
    if (any(grepl("nil",as.character(plot_levels)))) {
      plot_levels <- plot_levels[c(plot_levels[-which(plot_levels=="nil")], plot_levels[which(plot_levels=="nil")])]
    }
    # plot_levels <- factor(sort(nodes[nodes_module_entity$entity=="readout"]))
    plot_dynamics_df(elementary_states_dynamics, output_folder, plot_levels, font_size, plot_width, plot_height,line_thickness)
  }
  setwd("..")
  # probtraj_table_only_probs_TH
}

#######################################################################################

# create phenotype plots, extract table (df) of probability trajectories

create_plots_probtraj_table <- function(foldername, plot_width, plot_height, plotting_threshold, line_thickness, font_size) {
  
  output_filename <- paste(foldername, paste(gsub("^\\.","", foldername), "_probtraj_table.csv",sep=""), sep = "/")
  
  # if (length(unlist(strsplit(foldername,"/")))==1) {output_filename <- paste(foldername, gsub("\\/","", foldername), "_probtraj_table.csv", sep = "")}
  
  probtraj_table <- read.csv(output_filename, sep = "\t", stringsAsFactors = F,check.names=FALSE)
  if ( all(is.na(probtraj_table[,ncol(probtraj_table) ])) ) {
    probtraj_table <- probtraj_table[,-ncol(probtraj_table) ]
  }
  probtraj_table_only_probs <- probtraj_table[, c(1:4,  seq(5, ncol(probtraj_table), by=2) )]
  colnames(probtraj_table_only_probs) <- gsub("\\]","",gsub("Prob\\[","", colnames(probtraj_table_only_probs) ) )
  
  probtraj_table_only_probs <- probtraj_table_only_probs[,c(1,5:ncol(probtraj_table_only_probs))]
  
  columns_to_keep <- which( apply( probtraj_table_only_probs > plotting_threshold, 2, sum)>1  # has higher values than threshold
                            # & ( apply( probtraj_table_only_probs, 2, which.max) > 1 | # highest value is not first AND/OR last value not 0
                            #     probtraj_table_only_probs[nrow(probtraj_table_only_probs),] > plotting_threshold) 
                            ) 
  # probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>plotting_threshold
  
  probtraj_to_plot <- probtraj_table_only_probs[,columns_to_keep]
  
  # probtraj_table_only_probs[,c(1, which(probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>plotting_threshold)[
  # which(probtraj_table_only_probs[dim(probtraj_table_only_probs)[1],]>plotting_threshold)>4])]
  
  # if (ncol(probtraj_to_plot)>10) {
  #   prominent_states <- order(probtraj_to_plot[nrow(probtraj_to_plot),], decreasing = T)[1:10]
  #   prominent_states <- prominent_states[prominent_states!=1]
  #   probtraj_to_plot <- probtraj_to_plot[,c(1,prominent_states)]
  # }

  probtraj_table_only_nonzero_probs_df <- melt( probtraj_to_plot, id=c("Time"))
  
  # plot dynamics
  column_number <- round(sqrt(length(unique(probtraj_table_only_nonzero_probs_df$variable))))
  temporal <- ggplot(probtraj_table_only_nonzero_probs_df, aes(x=Time,y=value, color=variable)) + 
    geom_line(size=line_thickness) + ylab("Phenotype probability") +
    theme( legend.position = "right",legend.title = element_blank(), text=element_text(size = font_size),
           legend.text=element_text(size=font_size-2) ) +
    guides(fill = guide_legend(nrow = length(unique(probtraj_table_only_nonzero_probs_df$variable))))
  ggsave(filename = paste(foldername, "Phenotypes_probability_time_evolution.png", sep = "/"),
         plot=temporal, width = plot_width, height = plot_height)
  ####
  pie_data <- as.data.frame(t( probtraj_to_plot[ nrow(probtraj_to_plot), 2:ncol(probtraj_to_plot)]))
  colnames(pie_data)[1] <- "value"; pie_data$ord <- rownames(pie_data) # gsub("\\..","\\ & ",gsub("Prob.","",rownames(pie_data)))
  barplot_data <- pie_data[order(pie_data$value, decreasing = T),]
  barplot <- ggplot(barplot_data, aes(x = ord, y = value, label = value*100, fill=ord ) ) +
    geom_bar(stat = "identity") + geom_text(size = font_size/3, position = position_stack(vjust = 0.5)) + coord_flip() + 
    theme( axis.ticks=element_blank(), panel.grid=element_blank(),  # , axis.title= element_blank()
           legend.position = "none",legend.title = element_blank(), text=element_text(size=font_size) ) + ylab("Phenotype probability") + xlab("")
  ggsave(filename=paste(foldername, "Phenotypes_probability_barplot.png", sep = "/"), plot=barplot, width = plot_width, height = plot_height)
  
  # probtraj_table_only_probs[,probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>0]
  probtraj_to_plot
}

#######################################################################################

# this function splits "phenotypes' (simultaneously activated nodes in a sample trajectory) into their elementary states and sums these across simulations
from_phenotypes_to_state_probabilities <- function(df) {
  
  colindices <- lapply( unique(unlist(strsplit(colnames(df)[2:ncol(df)], "--"))), function(x) which(grepl(x, colnames(df))) )
  elementary_final_states <- cbind(df[,"Time"], data.frame(  matrix(0, nrow = nrow(df), ncol = length(colindices)) ))
  colnames(elementary_final_states) <- c("Time",unique(unlist(strsplit(colnames(df)[2:ncol(df)], "--"))))
  colnames(elementary_final_states) <- gsub(">","",gsub("<","",colnames(elementary_final_states)))
  
  for(counter in 1:length(colindices) ) {
    if (length(colindices[[counter]])>1) {
      elementary_final_states[,counter+1] <- rowSums(df[, colindices[[counter]]])
    } else {
      elementary_final_states[,counter+1] <- df[, colindices[[counter]]]
    }
  }
  if (any(grepl("nil",colnames(elementary_final_states)))){
    if ( all(elementary_final_states$nil[round(nrow(elementary_final_states)/2):nrow(elementary_final_states)]<10e-3) ) {
      elementary_final_states <- elementary_final_states[,-which(colnames(elementary_final_states)=="nil")]
    }
  }
  elementary_final_states
}

#######################################################################################

# create MABOSS models from table of BNET-format models

function_bnet_bnd_cfg_model_creator <- function(all_models_inputs_matrix,location_model_subfolders,path_bnet_to_maboss_py,variable_names,visible_nodes,
                                          filename_tag_list,nodes_to_set_initstate,initstate_values_list,fixed_genes,indices_fixed_genes,
                                          maboss_path,meta_sh_file_path){

lapply(1:dim(all_models_inputs_matrix)[2], function(x) {
  # generate bnet filename
  bnet_filename <- paste(location_model_subfolders,x,"/model",x,".bnet", sep="")
  
  # create bnet file, if it doesnt exist
  if (dir.exists(paste(location_model_subfolders,x,sep=""))==0){
    dir.create(paste(location_model_subfolders,x,sep=""))
  }

  if (file.exists(bnet_filename)==0){
    df <- data.frame(cbind(paste(variable_names,",",sep=""), all_models_inputs_matrix[,x]) )
    colnames(df) <- c("targets,", "factors")
    write_tsv( df, bnet_filename)
    # create bnd file, if it doesnt exist yet
    if (file.exists(gsub(".bnet",".bnd",bnet_filename))==0){
     # run python script bnet -> bnd
     system(paste("python ", path_bnet_to_maboss_py, rep("", sign(nchar(path_bnet_to_maboss_py))), 
                  "bnet_to_maboss.py ", bnet_filename, sep=""))
    }
  }

  # specify init states of CFG
  function_add_initial_states(gsub(".bnet",".cfg",bnet_filename), "", c(), "OFF", "")
  function_internal_visible_nodes(gsub(".bnet",".cfg",bnet_filename), visible_nodes)
  # CREATE cell-type + exp. condition specific CFG FILES
  function_cfg_file_creator(gsub(".bnet",".bnd", bnet_filename), gsub(".bnet",".cfg",bnet_filename), filename_tag_list,
                            nodes_to_set_initstate, initstate_values_list,
                            fixed_genes, indices_fixed_genes,"")
  
  # write(sh_command, sh_dir_file_name, app)
  
  # CREATE .sh files
  counter<-x
  sh_command_vars <- paste("cd ", maboss_path, "\n","source MaBoSS.env \n", 
                           "cd ", getwd(), "/",location_model_subfolders, counter, sep = "")
  
  cfg_pattern = glob2rx(paste("model", counter, "_*.cfg", sep = "" ))
  cfg_list <- list.files(paste(location_model_subfolders,counter,sep = ""), pattern=cfg_pattern, recursive = F)
  # generate list of commands
  sh_command <- paste(sh_command_vars,
                      paste0( paste("time MBSS_FormatTable.pl ", "model", counter,".bnd ", cfg_list, sep=""), collapse = "\n"), sep="\n")
  
  # write shell script with list of commands in folder
  sh_dir_file_name <- paste(location_model_subfolders,counter, "/run_simul.sh", sep="")
  write(sh_command, sh_dir_file_name)
  
  # RUN simulations
  # system(paste("cd", dir_name, "\n", "sh run_simul.sh") )
  
# write "meta" shell scripts that starts all other shell scripts
  if (!file.exists(meta_sh_file_path)){
    if (counter<dim(all_models_inputs)[2]){
      write( paste("sh krasmodel",counter,"/run_simul.sh & ", sep = ""), 
             paste(unlist(strsplit(location_model_subfolders,"/"))[[1]],"/run_all.sh",sep = ""), append=T)
    } else {
      write(paste("sh krasmodel",counter,"/run_simul.sh","\n","wait","\n","echo all complete", sep = ""), 
            "kras_model_topologies/run_all.sh", append=T)
    }
  }
})

}

#################################


create_probtraj_table <- function(foldername) {
  
  output_filename <- paste(foldername, paste(gsub("^\\.","", foldername), "_probtraj_table.csv",sep=""), sep = "/")
  
  # if (length(unlist(strsplit(foldername,"/")))==1) {output_filename <- paste(foldername, gsub("\\/","", foldername), "_probtraj_table.csv", sep = "")}
  
  probtraj_table <- read.csv(output_filename, sep = "\t", stringsAsFactors = F,check.names=FALSE)
  if ( all(is.na(probtraj_table[,ncol(probtraj_table) ])) ) {
    probtraj_table <- probtraj_table[,-ncol(probtraj_table) ]
  }
  probtraj_table_only_probs <- probtraj_table[, c(1:4,  seq(5, ncol(probtraj_table), by=2) )]
  colnames(probtraj_table_only_probs) <- gsub("\\]","",gsub("Prob\\[","", colnames(probtraj_table_only_probs) ) )
  
  probtraj_table_only_probs <- probtraj_table_only_probs[,c(1,5:ncol(probtraj_table_only_probs))]
  
  columns_to_keep <- which( apply( probtraj_table_only_probs > plotting_threshold, 2, sum)>1  # has higher values than threshold
                            # & ( apply( probtraj_table_only_probs, 2, which.max) > 1 | # highest value is not first AND/OR last value not 0
                            #     probtraj_table_only_probs[nrow(probtraj_table_only_probs),] > plotting_threshold) 
  ) 
  # probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>plotting_threshold
  
  probtraj_table_only_probs[,columns_to_keep]
  
  # probtraj_table_only_nonzero_probs_df <- melt( probtraj_to_plot, id=c("Time"))
  # probtraj_table_only_probs[,probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>0]
  # probtraj_to_plot
}

#################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########################

# very slow for large dataframes! investigate why...
function_extract_maxpoint <- function(vars, all_simulation_timecourses_df) {
  
  z=all_simulation_timecourses_df; z[is.na(z)]=0
  result <- data.frame(matrix(NA,length(unique(all_simulation_timecourses_df$condition_folder)), length(vars)))
  colnames(result) <- vars
  
  z <- z %>% gather(Description, value, cc_input:cell_death)
  
  for (variable in vars) {
    result_tmp <- z %>% group_by(condition_folder) %>% dplyr::filter(Description %in% variable) %>% dplyr::filter(value == max(value)) %>%
          group_by(condition_folder) %>% dplyr::filter(Time == min(Time))
    result[,variable] <- result_tmp$value
  }
  result
}

##########################
# convert boolnet model to SIF table

function_boolnet_to_SIF <- function(boolnet_table){
 if (sum(class(boolnet_table)!="data.frame")>0) {
    boolnet_table <- as.data.frame(boolnet_table)
  }
  # if there are parenth with ! before them, these have to be 
  for (counter in 1:nrow(boolnet_table)) {
    par_char <- boolnet_table$factors[counter]
    # "!(A | B) | !(C | D)" should become "!A | !B | !C | !D" (this is not exactly correct in terms of logic, but we just want to get influence graph)
    if (length(unlist(str_match_all(par_char, "!\\(.*?\\)")))>0){
      par_substr <- unlist(str_match_all(par_char, "!\\(.*?\\)"))
      for (substr_char in par_substr) {
        new_str <- gsub("!","",gsub(" ","",gsub("\\)","",gsub("\\(","",substr_char))))
        vars <- paste("!",unlist(strsplit(unlist(strsplit(new_str, "\\&")),"\\|")),sep=""); connectors <- unlist(str_match_all(new_str, "[^A-Za-z0-9$]"))
        new_str <- paste0(c(vars,connectors)[order(c(seq_along(vars), seq_along(connectors)))], collapse = "")
        par_char <- gsub(par_substr[par_substr %in% substr_char],new_str, par_char,fixed = T)
      }
      boolnet_table$factors[counter] <- par_char
      # print(counter)
    }
  }
  # remove spaces
  boolnet_table$factors <- gsub(" ","",boolnet_table$factors)
  z <- boolnet_table %>% mutate(factors = strsplit(factors, "&")) %>% unnest(factors) %>% mutate(factors = strsplit(factors, "\\|")) %>% unnest(factors)
  z$interaction <- 1; z$interaction[grepl("!",z$factors)] <- -1
  z$factors <- gsub("!","",z$factors); 
  # parentheses can be removed if there is no exclamation mark before
  z$factors <- gsub("\\)","",z$factors); z$factors <- gsub("\\(","",z$factors)
  z <- z[!duplicated(z),]; rownames(z)<-c()
  z[,c(2,3,1)]
}

##########################

# path <- paste(metafolder,"/",subfolder_stem, model_id,"/model", model_id,".bnet",sep="")
function_extract_bnet_convert_to_sif <- function(path) {
  boolnet_table <- as.data.frame(read_tsv(path))
  colnames(boolnet_table) <- c("targets", "factors")
  boolnet_table$targets <- gsub(",","",boolnet_table$targets)
  function_boolnet_to_SIF(boolnet_table)
}

##########################

function_model_table_to_redundant_sif <- function(all_models_inputs_matrix,variable_names,filepath,filename) {
  
  lapply(1:dim(all_models_inputs_matrix)[2], function(x) {
    df <- data.frame(cbind(paste(variable_names,",",sep=""), all_models_inputs_matrix[,x]) )
    colnames(df) <- c("targets", "factors"); rownames(df)<-c(); 
    df[,colnames(df)[1]] <- as.character(df[,colnames(df)[1]]); df[,colnames(df)[2]] <- as.character(df[,colnames(df)[2]])
    df$targets <- gsub(",","",df$targets)
    sif_table <- cbind(function_boolnet_to_SIF(df), model_ids[x]); colnames(sif_table)[length(colnames(sif_table))] <- "model_id"
    write_tsv(sif_table, paste(filepath,filename,sep=""), append=(!x==1)) # col_names = F,
  } )
}

##########################
##########################

# METAFUNCTIONS
# max points
function_maxpoints <- function(all_simulation_timecourses_df, all_simulation_endpoints, visible_nodes, max_variables,time_colname){
z <- all_simulation_timecourses_df; z[is.na(z)]=0
# max mk2
for (varname in max_variables) {
  # maximal chek1
  if (varname==max_variables[1]){
  result <- z %>% group_by(model_folder,condition_folder) %>% dplyr::filter(get(varname) == max(get(varname)) ) %>% 
    group_by(model_folder,condition_folder) %>% dplyr::filter(Time == min(Time))
  } else {
  result[,varname] <- (z %>% group_by(model_folder,condition_folder) %>% dplyr::filter(get(varname) == max(get(varname))) %>% 
                                  group_by(model_folder,condition_folder) %>% dplyr::filter(Time == min(Time)))[,varname]
  }
}

# maxpoints dataframe
rm(z); all_simulation_maxpoints <- all_simulation_endpoints; 
all_simulation_maxpoints[,max_variables] <- result[,max_variables]; all_simulation_maxpoints[,time_colname]<-result[,time_colname]
all_simulation_maxpoints <- all_simulation_maxpoints %>% gather(variable, value, visible_nodes)
all_simulation_maxpoints[is.na(all_simulation_maxpoints)]<-0
all_simulation_maxpoints
}

##########
# format exp values
function_format_expvalues <- function(exp_values,first_column_name, conditions_list){
colnames(exp_values) <- c()
if (is.character(exp_values[,1])) {
  rownames(exp_values) <- exp_values[,1]; exp_values <- exp_values[order(rownames(exp_values)),]
  exp_values <- cbind(conditions_list, as.data.frame(t(exp_values[,-1]))) # [match(1:length(conditions_list),conditions_order)][condition_id]
  colnames(exp_values)[1] <- first_column_name
    exp_values <- exp_values %>% gather(variable, value, colnames(exp_values)[!colnames(exp_values) %in% first_column_name])
}
exp_values
}

##########
# create empty DF to store simul results in

# column_names<-c("Time", "model_folder", "condition_folder", visible_nodes)
# all_simulation_timecourses_df <- function_create_simul_empty_df(all_cond_subfolders,rownumber,timestep_num, timevector, column_names)

function_create_simul_empty_df <- function(all_cond_subfolders,rownumber,timestep_num, timevector, 
                                           column_names, modelfoldername_replace, modelnumber_repl_regex){
rownumber <- length(all_cond_subfolders) 
# (length(list.dirs(metafolder,recursive=T))-1) - length(list.dirs(metafolder,recursive=F))
# length(list.dirs(".",recursive=T)) - (length(list.dirs(".",recursive=F))+1)

all_simulation_timecourses_df <- data.frame(matrix(NA, nrow=rownumber*timestep_num, ncol=length(column_names)))
colnames(all_simulation_timecourses_df) <- column_names
all_simulation_timecourses_df[,column_names[1]] <- rep(timevector,rownumber)
depth_model_folders <- unique(sapply(strsplit(all_cond_subfolders,"/"), length))-1
all_simulation_timecourses_df[,column_names[2]] <- unlist(lapply(sapply(strsplit(all_cond_subfolders,"/"),`[`,depth_model_folders), 
                                                            function(x) rep(x,timestep_num)))
# cut off string "krasmodel"
all_simulation_timecourses_df[,column_names[2]] <- as.numeric(gsub(modelfoldername_replace,"", 
                                                                   as.vector(as.matrix(all_simulation_timecourses_df[,column_names[2]]))))
# condition names
depth_cond_folders <- unique(sapply(strsplit(all_cond_subfolders,"/"), length))
all_simulation_timecourses_df[,column_names[3]] <- unlist(lapply(sapply(strsplit(all_cond_subfolders,"/"),`[`, depth_cond_folders),
                                                                function(x) rep(x,timestep_num)))
# cut off string "model"
all_simulation_timecourses_df[,column_names[3]] <- gsub(modelnumber_repl_regex,"",
                                                        as.vector(as.matrix(all_simulation_timecourses_df[,column_names[3]])) )
all_simulation_timecourses_df
}

#########

# function to extract all simuls

function_extract_simuls <- function(all_cond_subfolders, all_simulation_timecourses_df, variable_names,
                                    timestep_num, directory_depth, directory_index){
for (x in all_cond_subfolders) {
  # enter model folder
  setwd(apply(sapply(strsplit(x,"/"),`[`,directory_depth), directory_index, paste, collapse='/'))
  # extract from condition-specific subfolder
  probtraj_table_only_probs_TH <- create_probtraj_table(unlist(strsplit(x,"/"))[length(unlist(strsplit(x,"/")))])
  elementary_states <- from_phenotypes_to_state_probabilities(probtraj_table_only_probs_TH)
  row_indices <- seq( (which(all_cond_subfolders %in% x)-1)*timestep_num+1,
                      which(all_cond_subfolders %in% x)*timestep_num,1 )
  if (sum(colnames(elementary_states[,2:ncol(elementary_states)]) %in% colnames(all_simulation_timecourses_df))){
    elementary_states <- elementary_states[1:timestep_num, colnames(elementary_states) %in% variable_names]
    all_simulation_timecourses_df[row_indices, colnames(elementary_states)] <- elementary_states
  }
  print(x)
  setwd("../..")
}
  
  all_simulation_timecourses_df[,colnames(all_simulation_timecourses_df) %in% variable_names] <- 
    round(all_simulation_timecourses_df[,colnames(all_simulation_timecourses_df) %in% variable_names],3)
  all_simulation_timecourses_df
}

########


# function reorder allcond folders

function_reorder_allcond_folders <- function(subfolder_dirs, all_cond_subfolders, conditions_order, timestep_num){
  
# names of different conditions  
  conditions_list <- unique(gsub("model[0-9]+_","", unlist(lapply(
    sapply(strsplit(all_cond_subfolders,"/"),`[`,unique(sapply(strsplit(all_cond_subfolders,"/"), length))), 
    function(x) rep(x,timestep_num)))))
  
# re-order by MODEL folder (#2 should be before 11)
all_cond_subfolders <- all_cond_subfolders[order(as.numeric(gsub("[a-z]","",sapply(strsplit(all_cond_subfolders,"/"),`[`,3) )))]
# reorder condition names too
all_cond_subfolders <- paste(apply(sapply(strsplit(all_cond_subfolders,"/"),`[`,1:3), 2, paste, collapse='/'), 
                             paste(rep(paste0("model", 1:length(subfolder_dirs) ), each=length(conditions_list)),
                                   rep(conditions_list[match(1:length(conditions_list),conditions_order)],
                                       length(subfolder_dirs) ),sep="_"),sep="/")
# add dot before name
if (sum(substr(all_cond_subfolders,1,1)==".")<length(all_cond_subfolders)){
  all_cond_subfolders <- paste(".",all_cond_subfolders,sep="")
}
all_cond_subfolders
}

########

# DF for plotting results of 1 simul

# df_plot <- function_df_indiv_heatmap(all_simulation_maxpoints,vars_to_plot,model_id,exp_values)

function_df_indiv_heatmap <- function(all_simulation_maxpoints,vars_to_plot,model_id,exp_values){
df_plot <- all_simulation_maxpoints[all_simulation_maxpoints$variable %in% vars_to_plot &
                                      all_simulation_maxpoints$model_folder %in% model_id,]
df_plot$mse <- paste("MSE=", round(mse[df_plot$model_folder],3) ,sep="")
df_plot$model_folder <- paste("model #", df_plot$model_folder, sep = "")
#
# binding df with data
df_plot <- bind_rows(df_plot,exp_values); df_plot$mse[is.na(df_plot$mse)] <- "data"
df_plot$model_folder[is.na(df_plot$model_folder)] <- ""
df_plot
}

# create subfolder string
function_subfolder_string <- function(metafolder,subfolder_stem) {
  list.dirs(gsub(paste(getwd(),"/",sep=""),"",metafolder),recursive=F)[
  grepl(subfolder_stem, list.dirs(gsub(paste(getwd(),"/",sep=""),"",metafolder),recursive=F))]
}

# create sh with all commands
function_meta_sh_string <- function(subfolder_stem,subfolder_dirs,sh_file_name){
  meta_sh_string <- paste("sh ", subfolder_stem, seq(1,length(subfolder_dirs),1),sh_file_name, " &", sep = "");
  meta_sh_string[length(meta_sh_string)] <- gsub(" &","",meta_sh_string[length(meta_sh_string)])
  meta_sh_string
}

############

# crearing DF for indiv dynamics

function_indiv_model_dynamics <- function(all_simulation_timecourses_df,model_id,model_folder_colname,cond_folder_colname,conditions_list,vars_to_plot, time_colname){
  extracted_row_indices <- which(as.vector(as.matrix(all_simulation_timecourses_df[,model_folder_colname])) %in% model_id & 
                                   as.vector(as.matrix(all_simulation_timecourses_df[,cond_folder_colname])) %in% conditions_list )
  # conditions_list[match(1:length(conditions_list),conditions_order)][condition_id]
  # LOAD data (check variables)
  if (all(vars_to_plot %in% colnames(all_simulation_timecourses_df))) {
    indiv_model_dynamics <- all_simulation_timecourses_df[extracted_row_indices,
                                c(colnames(all_simulation_timecourses_df)[1:3],vars_to_plot)] %>% gather(variable, value,vars_to_plot)
    if (class(indiv_model_dynamics$value)!="numeric") {
      indiv_model_dynamics$value <- as.numeric(indiv_model_dynamics$value)
    }
    indiv_model_dynamics[as.vector(as.matrix(indiv_model_dynamics[,time_colname]))<time_limit,]
  } else {print("BAD VARIABLE NAMES!")}
}


######################

# select to plot value or error of variables

function_col_select <- function(df_plot, what_to_plot) {
  if (grepl("vals",what_to_plot)) { df_plot$value <- df_plot$value.x 
  if (grepl("round",what_to_plot)) { df_plot$value <- round(df_plot$value) }
  } else { df_plot$value <- abs(df_plot$value.x - df_plot$value.y)
  if (grepl("round",what_to_plot)) { df_plot$value <- round(df_plot$value) }
  }
  df_plot
}

############################

function_save_best_models <- function(all_models_inputs,smallest_n,mse,name_best5_sif,metafolder) {
  
  best_models_save <- cbind(all_models_inputs[order(mse)[1:smallest_n],], order(mse)[1:smallest_n], round(mse[order(mse)[1:smallest_n]],3))
  colnames(best_models_save)[(ncol(best_models_save)-1):ncol(best_models_save)] <- c('model_index', 'MSE')
  write_csv(best_models_save,paste(metafolder, "best5_models_logical.csv",sep="/"))
  # number of distinct rules per node
  # sapply(apply(all_models_inputs[order(mse)[1:smallest_n],], 2, unique),length)
  #
  # variable RULEs
  # all_models_inputs[order(mse)[1:smallest_n],
  #                  sapply(apply(all_models_inputs[order(mse)[1:smallest_n],], 2, unique),length)>1]
  
  # to visualize with caspots need boolnet model
  
  for (model_counter in 1:smallest_n) {
    model_id <- order(mse)[model_counter]
    path <- paste(metafolder,"/",subfolder_stem,model_id,"/model", model_id,".bnet", sep = "")
    best_model_sif <- function_extract_bnet_convert_to_sif(path) # 
    # SIF table with all models
    write_tsv(cbind(best_model_sif,model_id),paste(metafolder,"/",name_best5_sif,sep=""),append=(!model_counter==1))
  }
  
  # collapse identical rows, with model indices and # of models summarized 
  # (freq columns showing fraction of top 5 models where edge is present)
  best_model_sif_top5 <- read_tsv(paste(metafolder,"/",name_best5_sif,sep="")) %>% group_by(factors,interaction,targets) %>% 
    summarise_all(funs(paste(model_id,collapse=",")))
  # add frequencies
  best_model_sif_top5$freq <- sapply(strsplit(best_model_sif_top5$model_id,"\\,"),length)/length(
    unique(unlist(strsplit(best_model_sif_top5$model_id,"\\,"))))
  # SAVE
  write_tsv(best_model_sif_top5[,c(3:1,4:ncol(best_model_sif_top5))], paste(metafolder,"/best5_models_freq.sif",sep=""))
}

#########################################################

function_calc_mse <- function(all_simulation_maxpoints, exp_values){
  
  simul_fit_results <- all_simulation_maxpoints[all_simulation_maxpoints$variable %in% unique(exp_values$variable),] %>%
    arrange(model_folder, variable)
  
  # calculate difference
  if (sum( simul_fit_results$variable[simul_fit_results$model_folder==1] == exp_values$variable )==nrow(exp_values) 
      # sum(simul_fit_results$variable %in% exp_values$variable)==nrow(simul_fit_results) & 
      ) {
    if (class(simul_fit_results$value)!="numeric") {
      simul_fit_results$value <- as.numeric(simul_fit_results$value)
    }
    # mse <- 
    colSums(matrix( ( simul_fit_results$value - rep(exp_values$value, nrow(simul_fit_results)/nrow(exp_values) ) )^2, 
                   nrow(exp_values)))/nrow(exp_values)
    
  } else {
    print("mismatch in order of variables")
  }
}

##########################################################

# write bnet files from model table generated in MATLAB
# matlab_filepath <- "/data/users/mkoltai/research/models/KRAS_DNA_repair_model/matlab_ode/"
# model_ids <- c(20, 4, 12, 2, 18)
# all_models_inputs <- read.csv(paste(matlab_filepath,"best5_models_bnet.csv",sep = ""), stringsAsFactors = F) # as.data.frame()
# # create redundant SIF table
# function_model_table_to_redundant_sif(t(all_models_inputs), variable_names, matlab_filepath, "best5_models_redundant.sif")
# # create concise SIF table with frequencies
# best_model_sif_top5 <- read_tsv(paste(matlab_filepath,"/best5_models_redundant.sif",sep="")) %>% 
#   group_by(targets,interaction,factors) %>% summarise(model_ids = paste(model_id, collapse ="_"), freq=length(model_id)/smallest_n)
# # write CSV
# write_csv(best_model_sif_top5, paste(matlab_filepath,"/best5_models_frequencies.csv",sep=""))

########################################################################

# this is the initial pipeline for a single model version under different conditions

# # sample code to use functions
# 
# # nodes w specified initial states
# # define manually the tags to be added to the different model versions, as they are defined in "initstate_values_list". 
# filename_tag_list <- c("_wt","_DNAdam","_KRASmut")
# # initstate lists: the nodes whose initial states are defined
# nodes_to_set_initstate <- c("KRAS","DNA_damage")
# # the values that in different versions of the model these nodes assume
# initstate_values_list <- list(c(0,0), c(0,1), c(1,0));
# 
# # FIXED GENES
# input_cfg_file <- "kras_model_conceptual.cfg"; input_bnd_file <- gsub(".cfg",".bnd",input_cfg_file)
# # to what values these genes would be fixed: -1 is KO, +1 knock-in
# fixed_genes <- c("CHEK1","MAPKAPK2"); indices_fixed_genes <- list(c(0,-1),c(-1,0),c(-1,-1))
# 
# function_cfg_file_creator(input_cfg_file, filename_tag_list, nodes_to_set_initstate, initstate_values_list, fixed_genes, indices_fixed_genes)
# 
# ########################################################################
# # RUN SIMULATION
# 
# # SPECIFY path to MaBoSS.env and MBSS_FormatTable.pl (set it accordingly)
# dir_name <- gsub(".cfg","",input_cfg_file);
# sh_command_vars <- paste("cd ../.. \ncd MaBoSS-env-2.0/ \nsource MaBoSS.env \ncd .. \ncd", 
#                          paste(strsplit(getwd(), "/")[[1]][sapply(strsplit(getwd(), "/"), length)], dir_name, sep = "/") )
# 
# # generate list of commands
# sh_command <- paste(sh_command_vars, paste0(paste("time MBSS_FormatTable.pl", input_bnd_file, cfg_list ), collapse = "\n"), "cd ..", sep = "\n")
# # write shell script with list of commands in folder
# write(sh_command, paste(dir_name,"run_simul.sh",sep="/"))
# # RUN simulations
# system(paste("cd", dir_name, "\n", "sh run_simul.sh") )
# # system("rm run_simul.sh")
# 
# ########################################################################
# # PLOTTING
# # GET the folder that was changed in last x minutes (set the value accordingly)
# time_window <- 10; 
# # get the folder names with the data we want to plot
# output_folders <- system(paste("find ./*/* -type d -mmin -", time_window, " -print", sep = ""), intern=TRUE); output_folders <- gsub("/.*/","/",output_folders)
# # PLOTS
# # in each folder 4 plots are generated:
# # elementary_states_dynamics_barplot.png: dynamics of elementary states (nodes being ON)
# # elementary_states_dynamics_time_evolution.png: steady state of elementary states (nodes being ON)
# # Phenotypes_probability_barplot.png: dynamics of "phenotypes" (simultaneously activated nodes)
# # Phenotypes_probability_time_evolution.png: steady state of "phenotypes" (simultaneously activated nodes)
# #
# # parameters for plots. Plotting threshold is the minimal relative frequency of a single state (a node being ON) for elementary_states_*.png
# font_size <- 7; plot_width <- 6; plot_height <- 4; font_size<-12; plotting_threshold <- 0.001; line_thickness <- 1
# multiple_folder_plotter(dir_name, output_folders, plot_width, plot_height, plotting_threshold, line_thickness, font_size)
# 
# # plot and access data of an individual simulation (instead of all of one)
# # which folder?
# index_folder <- 1; plot_folder <- output_folders[index_folder] # paste(gsub(".bnd","",input_bnd_file),gsub("\\.","",output_folders[index_folder]),sep="")
# # enter directory of meta-folder
# setwd( gsub(".bnd","",input_bnd_file) )
# # plot phenotypes and load probability trajectories of selected simulation into a dataframe
# probtraj_table_only_probs_TH <- create_plots_probtraj_table(plot_folder, plot_width, plot_height, plotting_threshold,line_thickness,font_size)
# setwd("cd ..")
# 
# # access elementary states of phenotypes
# elementary_states <- from_phenotypes_to_state_probabilities(probtraj_table_only_probs_TH)
# 
# ########################################################################
# # create & simulate MaBoSS model from boolnet
# 
# # toymodel_cellfate.bnet
# 
# system("python bnet_to_maboss.py toymodel_cellfate.bnet")
# 
# # create range of models
# 
# # set up init conditions
# input_cfg_file <- "toymodel_cellfate.cfg"; input_bnd_file <- gsub(".cfg",".bnd",input_cfg_file); dir_name <- gsub(".cfg","",input_cfg_file);
# initstate_values_list <- list(0.1,0.5,0.7,1); nodes_to_set_initstate <- "A"
# filename_tag_list <- paste("A_initst", unlist(initstate_values_list), sep = "")
# 
# function_cfg_file_creator(input_cfg_file, filename_tag_list, nodes_to_set_initstate, initstate_values_list, list(), list())
# 
# ########################################################################
# # SIMULATE
# 
# go_maboss_path <- "cd /bioinfo/users/mkoltai/research/models/MaBoSS-env-2.0/"
# input_bnd_file <- gsub(".cfg",".bnd",input_cfg_file); dir_name <- gsub(".cfg","",input_cfg_file);
# simulator_function(input_cfg_file, input_bnd_file, dir_name, go_maboss_path)
# 
# ########################################################################
# # PLOT, ANALYSE
# time_window <- 5; 
# # get the folder names with the data we want to plot
# output_folders <- system(paste("find ./*/* -type d -mmin -", time_window, " -print", sep=""), intern=TRUE); output_folders <- gsub("/.*/","/",output_folders)
# font_size <- 7; plot_width <- 6; plot_height <- 4; font_size <- 12; plotting_threshold <- 0.01; line_thickness <- 1
# multiple_folder_plotter(dir_name, output_folders, plot_width, plot_height, plotting_threshold, line_thickness, font_size)

##########################################################

# version of the model with cell death deactivating repair and cdc25b, that yielded the best fit

# variable_names=c("cc","kras", "dna_dam", "chek1", "mk2", "atm_atr", "hr",
#                  "cdc25b", "g2m_trans", "cell_death") # 'phosphatase','proliferation'
# 
# # all inputs
# list_all_inputs <- list(
#   cc_inputs="cc",
#   kras_inputs=c('kras'),
#   dna_dam_inputs=c('(dna_dam | kras) & !hr'), #  & !phosphatase
#   chek1_inputs=c('atm_atr'),
#   mk2_inputs=c('atm_atr | kras', 'atm_atr & kras'),
#   atm_atr_inputs=c('(dna_dam)','(dna_dam | chek1)'),
#   hr_inputs=c('(atm_atr | hr) & !cell_death'),
#   cdc25b_inputs = c('(cc | kras) & (!chek1 & !mk2)', '(cc | kras) & (!chek1 & !mk2) & !cell_death'), # !cell_death  & !phosphatase'
#   g2m_trans_inputs=c('g2m_trans | cdc25b'),
#   cell_death_inputs=c('cell_death | (dna_dam & g2m_trans)') )


########################
# extract 1 simulation
# setwd("kras_model_topologies/krasmodel120")
# probtraj_table_only_probs_TH <- create_probtraj_table("model120_KRASmut_fixed_mk2off")
# elementary_states <- from_phenotypes_to_state_probabilities(probtraj_table_only_probs_TH)
# setwd("../..")
######

# all inputs
# c("cc","kras", "dna_dam", "chek1", "mk2", "atm_atr", "hr", "cdc25b", "g2m_trans", "cell_death") # , 'phosphatase','proliferation'
# list_all_rules <- list(
#         cc_inputs="cc",
#         kras_inputs=c('kras'),
#         dna_dam_inputs=c('(dna_dam | kras) & !hr'), #  & !phosphatase
#         chek1_inputs=c('atm_atr'),
#         mk2_inputs=c('atm_atr | kras', 'atm_atr & kras'),
#         atm_atr_inputs=c('(dna_dam)','(dna_dam | chek1)'),
#         hr_inputs=c('(atm_atr | hr) & !cell_death'),
#         cdc25b_inputs = c('(cc | kras) & (!chek1 & !mk2)', '(cc | kras) & (!chek1 & !mk2) & !cell_death'), # !cell_death  & !phosphatase'
#         g2m_trans_inputs=c('g2m_trans | cdc25b'),
#         cell_death_inputs=c('cell_death | (dna_dam & g2m_trans)') )
# phosphatase_inputs=c('cell_death'),
# proliferation='g2m_trans & !cell_death'
#
# best performing model currently is
# kras_model_topologies/results/krasmodel_repair_one_variable/dnadam_atmatr_chk1mk2_hr_celldeathdeactivation

# best model
# best_models_logical <- as.data.frame(read_csv(
# "~/research/models/KRAS_DNA_repair_model/maboss_models/kras_model_topologies/results/krasmodel_repair_one_variable/dnadam_atmatr_chk1mk2_hr_celldeathdeactivation/best5_models_logical.csv"))


###########################################

function_convert_paramscan_to_table <- function(filename){
  parscan = unlist(strsplit(read_file(filename),"\n"))
  parscan=parscan[which(grepl("Parameter set", parscan))[1]:length(parscan)]; parscan = parscan[!parscan %in% ""]
  
  paramsets_df = as.data.frame(t(as.data.frame(strsplit(gsub(" ","",gsub("\\$","",gsub("Parameter set : ","",parscan[grepl("Parameter set",parscan)]))), ","))))
  colnames(paramsets_df)=unlist(lapply(1:ncol(paramsets_df), function(n) {unique(gsub("=[0-9]","", paramsets_df[,n] ))}))
  rownames(paramsets_df)=c(); paramsets_df[]=lapply(paramsets_df, function(x) as.numeric(gsub("[^0-9]","",x)))
  
  parscan_df = as.data.frame(t(as.data.frame(sapply(parscan[grepl("#",parscan)], function(x) {strsplit(x,"\t")})))); rownames(parscan_df)=c()
  parscan_df$paramset_id = NA; parscan_df$paramset_id[parscan_df$V2 %in% "Value"] <- 1:sum(parscan_df$V2 %in% "Value")
  parscan_df = parscan_df %>% fill(paramset_id)
  # fill in parameter values
  parscan_df[,colnames(paramsets_df)]<-paramsets_df[parscan_df$paramset_id,]
  # remove redundant rows
  colnames(parscan_df)[2:(which(colnames(parscan_df) %in% "paramset_id")-1)]=as.vector(as.matrix(parscan_df[1,2:(which(colnames(parscan_df) %in% "paramset_id")-1)]))
  parscan_df = parscan_df[!(parscan_df[,2] %in% "Value"),]
  # to get average, multiply values by frequency
  value_indices=3:(which(colnames(parscan_df) %in% "paramset_id")-1)
  parscan_df[,value_indices] = as.data.frame( apply(parscan_df[,value_indices], 2, as.numeric)*as.numeric(as.vector(as.matrix(parscan_df$Value))) )
  par_cols=colnames(parscan_df)[(which(colnames(parscan_df) %in% "paramset_id")+1):ncol(parscan_df)]
  parscan_df[,3:ncol(parscan_df)] %>% group_by(paramset_id) %>% summarise_all(funs(sum))
}


################################################

function_meanvalues_paramscan <- function(paramscan_path,filename_parscan){
  pscan_results=read_tsv(paste(paramscan_path,filename_parscan,sep=""))
 # columns with variables
  value_indices = (which(colnames(pscan_results) %in% "Value")+1):(min(which(grepl("\\$",colnames(pscan_results))))-1)
 # weight by frequency
 pscan_results[,value_indices]=as.data.frame(apply(pscan_results[,value_indices],2,as.numeric)*as.numeric(as.vector(as.matrix(pscan_results$Value))) )
 # sum for each param set
  states_sum_df=pscan_results %>% group_by(`ParameterSet Id`) %>% summarise_at(.vars=names(.)[value_indices],.funs=sum)
  pscan_results_means=left_join(states_sum_df, unique(pscan_results[,colnames(pscan_results) %in% "ParameterSet Id" | grepl("\\$",colnames(pscan_results))]) )
 # pscan_results_means=pscan_results_means[,!colnames(pscan_results_means) %in% c("Steady state","Value")]
 # sort columns (parameters on the left)
  pscan_results_means[,sort(colnames(pscan_results_means))]
}

################################################

function_create_conditions_df_json <- function(nodes_to_set_initstate, initstate_values_list,fixed_genes,indices_fixed_genes){
  counter=0
  for (initstate_counter in 1:length(initstate_values_list)) {
    for (fixedparam_counter in 1:length(indices_fixed_genes)) {
      counter=counter+1
      # fixed params
      if (sum(unlist(indices_fixed_genes[fixedparam_counter])!=0)>0){
        m_fixed_params=matrix(NA, nrow=length(unlist(indices_fixed_genes[fixedparam_counter])!=0), ncol=2)
        # KOs
        if (sum(unlist(indices_fixed_genes[fixedparam_counter])<0)>0){
          m_fixed_params[unlist(indices_fixed_genes[fixedparam_counter])<0,1]=paste(rep("$u_",sum(unlist(indices_fixed_genes[fixedparam_counter])<0)), 
                                                                                    fixed_genes[unlist(indices_fixed_genes[fixedparam_counter])<0], sep="")
          m_fixed_params[unlist(indices_fixed_genes[fixedparam_counter])<0,2]=rep(0,sum(unlist(indices_fixed_genes[fixedparam_counter])<0))}
        # Kins
        if (sum(unlist(indices_fixed_genes[fixedparam_counter])>0)>0){
          m_fixed_params[unlist(indices_fixed_genes[fixedparam_counter])<0,1]=paste(rep("$d_",sum(unlist(indices_fixed_genes[fixedparam_counter])<0)), 
                                                                                    fixed_genes[unlist(indices_fixed_genes[fixedparam_counter])<0], sep="")
          m_fixed_params[unlist(indices_fixed_genes[fixedparam_counter])<0,2]=rep(0,sum(unlist(indices_fixed_genes[fixedparam_counter])<0))}
      } else { m_fixed_params=matrix(NA,nrow=1,ncol=2)   }
      
      # fixedgene matrix: nodes
      if (sum(unlist(indices_fixed_genes[fixedparam_counter])!=0)>0){
        vals=sign(unlist(indices_fixed_genes[fixedparam_counter])+1)[unlist(indices_fixed_genes[fixedparam_counter])!=0]
        fixedgene_matrix_nodes=cbind(fixed_genes[unlist(indices_fixed_genes[fixedparam_counter])!=0], vals)
      } else { fixedgene_matrix_nodes=matrix(NA,nrow=1,ncol=2) }
      
      # all conditions (for 1 cell line)
      df_condition=as.data.frame(rbind( cbind(nodes_to_set_initstate, initstate_values_list[[initstate_counter]]), fixedgene_matrix_nodes, m_fixed_params))
      colnames(df_condition)=c("name","value"); df_condition=df_condition[apply(is.na(df_condition),1,sum)==0,]
      df_condition$cond_id=counter;
      
      if (initstate_counter*fixedparam_counter==1) {
        df_all_conditions=df_condition
      } else {
        df_all_conditions=rbind(df_all_conditions,df_condition)
      }
      
    }
  }
  rownames(df_all_conditions)<-c(); 
  if (!is.numeric(df_all_conditions$value)){
    df_all_conditions$value=as.numeric(as.character(df_all_conditions$value))}
  if (!is.numeric(df_all_conditions$cond_id)){
    df_all_conditions$cond_id=as.numeric(as.character(df_all_conditions$cond_id))}
  df_all_conditions
  
}

###############################

function_extract_csv_param_sampling <- function(filelist,id_vars) {
  # read in files to one DF
  for (filelist_variable in filelist) {
    results_untidy_df=function_meanvalues_paramscan(paramscan_path,filelist_variable)
    results_untidy_df$cond_id=as.numeric(str_replace(strsplit(filelist_variable,"cond")[[1]][2],".csv","")) 
    if (filelist_variable==filelist[1]){
      pscan_results_means_all_conditions=results_untidy_df
    } else {
      pscan_results_means_all_conditions=rbind(pscan_results_means_all_conditions,results_untidy_df)
    }
  }
  # pscan_results_means_all_conditions
  melt(pscan_results_means_all_conditions,id.vars=id_vars)
}

###############

function_create_errors_df_parsampling <- function(pscan_results_tidy, exp_values, results_untidy_df, colnames_exp_values,name_paramset_id){
  
  if (class(pscan_results_tidy$variable)!="character") {
    pscan_results_tidy$variable=as.character(pscan_results_tidy$variable)
  }
  
  if (sum(colnames(pscan_results_tidy) %in% "condition_folder")==0){
    pscan_results_tidy=left_join(pscan_results_tidy, exp_values, by=colnames_exp_values)
    colnames(pscan_results_tidy) = str_replace(colnames(pscan_results_tidy),".x","_sim")
    colnames(pscan_results_tidy) = str_replace(colnames(pscan_results_tidy),".y","_exp")
  }
  
  # CALCULATE ERRORS
  pscan_results_tidy$abs_dev=abs(pscan_results_tidy$value_sim-pscan_results_tidy$value_exp)
  pscan_results_tidy$sq_dev=(pscan_results_tidy$value_sim - pscan_results_tidy$value_exp)^2
  # exp 
  exp_values_cond=exp_values[exp_values$cond_id %in% unique(pscan_results_tidy$cond_id),]
  
  mean_errors_paramset = pscan_results_tidy[pscan_results_tidy$variable %in% exp_values_cond$variable,] %>% 
    group_by(`ParameterSet Id`) %>% summarise_at(.vars = c("abs_dev","sq_dev"),.funs = mean)
  mean_errors_vars =     pscan_results_tidy[pscan_results_tidy$variable %in% exp_values_cond$variable,] %>% 
    group_by(variable) %>% summarise_at(.vars = c("abs_dev","sq_dev"),.funs = mean)
  
  if (sum(grepl("\\$",colnames(mean_errors_paramset)))==0){
    mean_errors_paramset = join(mean_errors_paramset, 
      unique(results_untidy_df[,colnames(results_untidy_df) %in% name_paramset_id | grepl("\\$",colnames(results_untidy_df))]),
      by=name_paramset_id)
    
    # rank by error
    mean_errors_paramset$rank_abs_dev=match(1:dim(mean_errors_paramset)[1],order(mean_errors_paramset$abs_dev))
    mean_errors_paramset$rank_sq_dev=match(1:dim(mean_errors_paramset)[1],order(mean_errors_paramset$sq_dev))
    
  }
  
  mean_errors_paramset
  
}


##############################

function_create_tidy_df_errors_parsampling <- function(mean_errors_paramset, errors, error_type_not_displayed){
  
  # met to tidy df
  mean_errors_paramset_tidy=melt(mean_errors_paramset,id.vars=joining_cols) 
  # REORDER factors
  mean_errors_paramset_tidy$variable=factor(mean_errors_paramset_tidy$variable,levels=sort(levels(mean_errors_paramset_tidy$variable)))
  
  ###################
  # PLOT of deviations per param set
  
  # errors=c("abs_dev", "sq_dev"); error_type_not_displayed=errors[1]; 
  error_disp=errors[!errors %in% error_type_not_displayed]
  # create a column where errors and also paramvals are rescaled to highlight relative differences
  mean_errors_paramset_tidy$value_mod=mean_errors_paramset_tidy$value
  mean_errors_paramset_tidy$value_mod[mean_errors_paramset_tidy$variable %in% error_disp] = 
    rescale(mean_errors_paramset_tidy$value[mean_errors_paramset_tidy$variable %in% error_disp])
  mean_errors_paramset_tidy$value_mod[grepl("\\$",mean_errors_paramset_tidy$variable)]=
    rescale(log10(mean_errors_paramset_tidy$value[grepl("\\$",mean_errors_paramset_tidy$variable)]))*(-1) # *maxval
  
  mean_errors_paramset_tidy
  
}

#####################################
# create JSON table

function_create_conditions_objectives_df = function(df_all_conditions,max_variables,exp_values,col_obj_states) {
  
  df_all_conditions$value=as.numeric(as.character(df_all_conditions$value))
  list_all_conditions = split(df_all_conditions[,colnames(df_all_conditions) %in% col_df_all_conditions], f = df_all_conditions$cond_id)
  # OBJECTIVES per cell line (from experimental data): 1st columns vars, 2nd value
  
  z=exp_values; colnames(z)[colnames(z) %in% "variable"]="name"; 
  z$type=NA; z$type[z$name %in% max_variables]="max"; z$type[!z$name %in% max_variables]="end"
  list_obj_states = split(z[,colnames(z) %in% col_obj_states], f=z$cond_id)
  
  # create dataframe whose columns are lists of 1) conditions 2) objectives
  num_cell_lines=length(list_obj_states)
  df_conditions_objectives=data.frame(matrix(NA,nrow=num_cell_lines,ncol=2))
  colnames(df_conditions_objectives)=colnames_df_conditions_objectives
  # fill in all conditions list
  df_conditions_objectives$conditions=list_all_conditions
  # fill in all objectives list
  df_conditions_objectives$objectives=list_obj_states
  df_conditions_objectives
}

#####################################

function_create_optim_param_df <- function(params_optim, init_vals, paramval_range, colnames_optim_param_df, digits_val) {
  optim_param_df=as.data.frame(matrix(NA,nrow=length(params_optim),ncol=length(colnames_optim_param_df))) 
  # json_df_settings$optimization_parameters
  colnames(optim_param_df)=colnames_optim_param_df; 
  optim_param_df$name=params_optim; optim_param_df$initial=init_vals; 
  optim_param_df$min=paramval_range[1]; optim_param_df$max=paramval_range[2]; optim_param_df$digits=digits_val
  optim_param_df
}

#####################################

function_sampling_results_cut <- function(sampling_results,max_variables, nodes_data_constr,id_cols,cond_col_name){
  # SUBSET for columns
  sampling_results_cut = sampling_results[,colnames(sampling_results)[grepl("\\$",colnames(sampling_results)) | 
                                                                        colnames(sampling_results) %in% c(id_cols, paste("max(",max_variables,")",sep = ""),
                                                                                                          nodes_data_constr[!nodes_data_constr %in% max_variables])]]
  colnames(sampling_results_cut) = gsub(id_cols[1],cond_col_name,gsub("\\)","",gsub("max\\(","",colnames(sampling_results_cut))))
  sampling_results_cut
}

#####################################

function_my_lines <- function(line_word,y_axis,plot_df) {
  if (grepl(line_word,y_axis)) {adj_val=0.5} else {adj_val=-0.5}; xval=length(unique(plot_df$variable[grepl("\\$",plot_df$variable)]))+0.5
  my_lines=data.frame(x=xval, y=adj_val, xend=xval, yend=length(unique(plot_df$y_axis))+adj_val)
  my_lines
}

#####################################

function_create_plot_df_parsampling_errors = function(mean_errors_paramset_tidy,cols_plot_df,first_col_name){
  plot_df=mean_errors_paramset_tidy[!mean_errors_paramset_tidy$variable %in% error_type_not_displayed, 
                                    colnames(mean_errors_paramset_tidy) %in% cols_plot_df]; colnames(plot_df)[1]=first_col_name
                                    plot_df
}

######################################

function_mse_plot <- function(mse, n_to_show,smallest_n,title_mse,xlab_var) {
  plot(order(mse)[1:n_to_show], mse[order(mse)[1:n_to_show]], xlab=xlab_var, ylab="",
       ylim=c(min(mse)*0.9,1.1*max(mse)), cex=1.5, cex.lab=1.5, cex.axis=1.2, xaxt="n" )
  axis(1,at=1:n_to_show,labels=as.character(1:n_to_show))
  title(ylab=title_mse, line=2.5, cex.lab=1.5, family="Calibri Light")
  points(order(mse)[1:smallest_n], mse[order(mse)[1:smallest_n]], pch=21, bg="red", col="red", cex=2.5)
  lapply(order(mse)[1:smallest_n], function(x) { text(x+0.6,mse[x]*1.05, as.character(x),cex=1.3) })  
}

############################################################
# sampling for a model w different cell lines (fixed param values and init conditions) 
# by running the sample separately for these conditions
# simulations launched by sh file "~/research/models/KRAS_DNA_repair_model/MaBoSS-Sampling/sampler_separate_runs.sh"
# containing commands: "./sampling -c *.cfg -p '$d_hr, $u_dna_dam, $u_cdc25b' -v '0.1, 1, 10' *.bnd -o *.csv"
# resulting simul files are in "~/research/models/KRAS_DNA_repair_model/maboss_models/param_scan/"
# WARNING: these simul files contain the last time points for all nodes
# filename_parscan="paramscan_all_conditions_endvals"
# paramscan_path="~/research/models/KRAS_DNA_repair_model/maboss_models/param_sampling_output/"
# # get file names, order file names by condition number
# filelist=list.files(paramscan_path,pattern="cond"); filelist=filelist[grepl(".csv",filelist)]
# filelist=filelist[order(as.numeric(str_replace(sapply(strsplit(filelist,"cond"),'[[',2),".csv","")))]
# # READ IN all CSVs
# id_vars=c("ParameterSet Id","cond_id")
# pscan_results_means_1cond=function_meanvalues_paramscan(paramscan_path,filelist[1]) # this is needed in another function
# pscan_results_means_all_conditions_tidy=function_extract_csv_param_sampling(filelist, id_vars)
# #### EXP vs SIMUL
# # LOAD exp_values (data)
# # exp_values <- as.data.frame(read_csv("../matlab_ode/exp_values_heatmap_rownames.csv",col_names = F))
# # # conditions_list = c("wt","wt_fixed_chek1off","wt_fixed_mk2off","wt_fixed_chek1off_mk2off",
# # #                     "DNAdam","DNAdam_fixed_chek1off","DNAdam_fixed_mk2off","DNAdam_fixed_chek1off_mk2off",
# # #                     "KRASmut","KRASmut_fixed_chek1off","KRASmut_fixed_mk2off","KRASmut_fixed_chek1off_mk2off")
# # exp_values <- function_format_expvalues(exp_values,"condition_folder",conditions_list)
# cond_ids=1:length(unique(exp_values$condition_folder))
# exp_values$cond_id=rep(cond_ids, nrow(exp_values)/length(cond_ids))
# 
# ### JOIN exp_values with SIMUL VALS
# colnames_exp_values = c("variable","cond_id"); name_paramset_id="ParameterSet Id"; joining_cols=c("ParameterSet Id","rank_abs_dev","rank_sq_dev")
# errors=c("abs_dev", "sq_dev"); error_type_not_displayed=errors[1]; 
# # DATAFRAME with errors
# mean_errors_paramset = function_create_errors_df_parsampling(pscan_results_means_all_conditions_tidy, exp_values, 
#                                                         pscan_results_means_1cond, colnames_exp_values,name_paramset_id)
# # TIDY DATAFRAME with errors for plotting
# mean_errors_paramset_tidy=function_create_tidy_df_errors_parsampling(mean_errors_paramset,errors,error_type_not_displayed)
# ###############
# # PLOT paramters
# #
# # to SAVE
# height_width_vals=c(10,8); paramscan_plot_name=filename_parscan # paste(gsub(".csv","",filename_parscan),"_all_cond",sep="")
# postscript(paste(paramscan_path, paramscan_plot_name,"u_dnadam_u_cdc25_d_hr_MEAN_ERROR.eps",sep=""), 
#            height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# ###############
# # plot by rank or by paramset index
# y_axis="rank_sq_dev"; cols_plot_df=c(y_axis,"variable","value","value_mod"); first_col_name="y_axis"
# plot_df=function_create_plot_df_parsampling_errors(mean_errors_paramset_tidy,c(y_axis,"variable","value","value_mod"),"y_axis")
# # demarc line
# line_word="rank"; my_lines=function_my_lines(line_word,y_axis,plot_df)
# ggplot(plot_df,aes(x=variable,y=y_axis,fill=value_mod)) + geom_tile(color="black",size=0.02) + geom_text(aes(label=round(value,2)), size=4) + 
#   scale_fill_gradient2(low="green",mid="white",high="red") + # scale_fill_gradient(low="white",high="red",limits=c(0,1)) + 
#   geom_segment(data=my_lines,aes(x,y,xend=xend, yend=yend),size=1.5,inherit.aes=F,color="blue") + ggtitle(paramscan_plot_name) + guides(fill=FALSE)
# # theme
# dev.off()

#####################

# theme1 <- theme(strip.background = element_blank(),
#                 axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
#                 strip.text.x=element_text(size=strip_text_x_size),strip.text.y=element_blank(), # strip.text.x=element_text(size=6)
#                 # element_blank()
#                 axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=axis_text_y_size), # 
#                 legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=legend_text_size), # element_blank()
#                 legend.box.background=element_rect(), legend.box.margin=margin(1,3,3,3), 
#                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())


######################################

function_boolnet_to_SIF <- function(boolnet_table){
  if (sum(class(boolnet_table)!="data.frame")>0) {
    boolnet_table <- as.data.frame(boolnet_table)
  }
  # if there are parenthesis with ! before them, these have to be 
  for (counter in 1:nrow(boolnet_table)) {
    par_char <- boolnet_table$factors[counter]
    # "!(A | B) | !(C | D)" should become "!A | !B | !C | !D" (this is not exactly correct in terms of logic, but we just want to get influence graph)
    if (length(unlist(str_match_all(par_char, "!\\(.*?\\)")))>0){
      par_substr <- unlist(str_match_all(par_char, "!\\(.*?\\)"))
      for (substr_char in par_substr) {
        new_str <- gsub("!","",gsub(" ","",gsub("\\)","",gsub("\\(","",substr_char))))
        vars <- paste("!",unlist(strsplit(unlist(strsplit(new_str, "\\&")),"\\|")),sep=""); connectors <- unlist(str_match_all(new_str, "[^A-Za-z0-9$]"))
        new_str <- paste0(c(vars,connectors)[order(c(seq_along(vars), seq_along(connectors)))], collapse = "")
        par_char <- gsub(par_substr[par_substr %in% substr_char],new_str, par_char,fixed = T)
      }
      boolnet_table$factors[counter] <- par_char
      # print(counter)
    }
  }
  # remove spaces
  boolnet_table$factors <- gsub(" ","",boolnet_table$factors)
  z <- boolnet_table %>% mutate(factors = strsplit(factors, "&")) %>% unnest(factors) %>% mutate(factors = strsplit(factors, "\\|")) %>% unnest(factors)
  z$interaction <- 1; z$interaction[grepl("!",z$factors)] <- -1
  z$factors <- gsub("!","",z$factors); 
  # parentheses can be removed if there is no exclamation mark before
  z$factors <- gsub("\\)","",z$factors); z$factors <- gsub("\\(","",z$factors)
  z <- z[!duplicated(z),]; rownames(z)<-c()
  z[,c(2,3,1)]
}

# provide path
function_extract_bnet_convert_to_sif <- function(path) {
  boolnet_table <- as.data.frame(read_csv(path))
  colnames(boolnet_table) <- c("targets", "factors")
  boolnet_table$targets <- gsub(",","",boolnet_table$targets)
  function_boolnet_to_SIF(boolnet_table)
}
