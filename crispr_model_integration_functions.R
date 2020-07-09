# functions used for crispri model integration

function_collapse_df_column_to_vector <- function(df) { as.vector(as.matrix(df)) }

#####################

function_find_row_interaction <- function(x,y,dataframe,source_column_number,target_column_number) {
  if (sum(class(dataframe)!="data.frame")>0){
    dataframe <- as.data.frame(dataframe)
  }
  # column indices for the pypath table
  if(missing(target_column_number)){
    target_column_number <- 3
  }
  if(missing(source_column_number)){
    source_column_number<-1; 
  }
  c(which( dataframe[,1] %in% x & dataframe[,3] %in% y ), which( dataframe[,1] %in% y & dataframe[,3] %in% x))
}

#####################

# two step paths

function_table_connecting_first_neighbors <- function(source,target,df_large_network,df_small_network,
                                                      source_column_number,target_column_number){
  # source and target nodes of small network  
  # first df of target neighbors
  df1 <- df_large_network[function_collapse_df_column_to_vector(df_large_network[,source_column_number]) %in% source,
                          source_column_number:target_column_number] 
  colnames(df1) <- c("source_hgnc","interaction_directed_signed1","intermed")
  df1 <- df1[!duplicated(df1) & !(df1[,source_column_number] == df1[,target_column_number]),]
  # second df of source neighbors
  df2 <- df_large_network[function_collapse_df_column_to_vector(df_large_network[,target_column_number]) %in% target,
                          source_column_number:target_column_number]
  colnames(df2) <- c("intermed","interaction_directed_signed2","target_hgnc")
  df2 <- df2[!duplicated(df2) & !(df2[,source_column_number] == df2[,target_column_number]),]
  # all paths source -> x -> target
  df_connecting_first_neighbors <- join(df1,df2,by="intermed",type="inner") # merge(df1,df2,by="intermed")
  df <- join(df_small_network, df_connecting_first_neighbors, by=c("source_hgnc", "target_hgnc"))
  df[apply(apply(df,2,is.na ),1,sum)<3,]
}


#######################

# two step paths

function_table_oneintermediate_paths <- function(source,target,df_large_network,source_column_number,target_column_number, colnames_table){
  
  df1 <- df_large_network[function_collapse_df_column_to_vector(df_large_network[,source_column_number]) %in% source,
                          source_column_number:target_column_number] 
  colnames(df1) <- c("source_hgnc","interaction_directed_signed1","intermed")
  df1 <- df1[!duplicated(df1) & !(df1[,source_column_number] == df1[,target_column_number]),]
  # second df of source neighbors
  df2 <- df_large_network[function_collapse_df_column_to_vector(df_large_network[,target_column_number]) %in% target,
                          source_column_number:target_column_number]
  colnames(df2) <- c("intermed","interaction_directed_signed2","target_hgnc")
  df2 <- df2[!duplicated(df2) & !(df2[,source_column_number] == df2[,target_column_number]),]
  # all paths source -> x -> target
  df_output <- join(df1,df2,by="intermed",type="inner")
  colnames(df_output) <- rev(rev(colnames_table)[1:length(colnames(df_output))])
  df_output
}

#####################

# gives three step (two intermediate) paths 

function_two_intermediate_nodes_path <- function(source, target, dataframe, source_column_number, target_column_number) {

  # if input table not dataframe (but tibble, typically)
  if (sum(class(dataframe)!="data.frame")>0){
    dataframe <- as.data.frame(dataframe)
  }

  # column indices for the pypath table
  if(missing(target_column_number)){
    target_column_number<-3
  }
  if(missing(source_column_number)){
    source_column_number<-1; 
  }

  source_first_neighbor <- dataframe[dataframe[,source_column_number] %in% source, target_column_number]
  target_first_neighbor <- dataframe[ dataframe[,target_column_number] %in% target, source_column_number]
  # source's 1st neighbour -> target's first neighbour
  paths_indices_x1_x2 <- which( dataframe[,source_column_number] %in% source_first_neighbor & dataframe[,target_column_number] %in% target_first_neighbor)
  
  c( # source -> x1
    which(dataframe[,target_column_number] %in% dataframe[paths_indices_x1_x2, source_column_number] & 
            dataframe[,source_column_number] %in% source),
    # x1 -> x2
    paths_indices_x1_x2,
    # x2 -> target
    which(dataframe[,source_column_number] %in% dataframe[paths_indices_x1_x2, target_column_number] & 
            dataframe[,target_column_number] %in% target) )
}

#####################
#####################

# function to find path with THREE intermediate nodes, extract central node
function_three_intermediate_nodes_central_node <- function(source,target,dataframe,source_column_number, target_column_number) {
  # column indices for the pypath table
  if(missing(target_column_number)){
    target_column_number <- 3
  }
  if(missing(source_column_number)){
    source_column_number<-1; 
  }
  
  intersect(
    # source of sources of target node
    dataframe[dataframe[,target_column_number] %in% dataframe[ dataframe[,target_column_number] %in% target, source_column_number], 
              source_column_number],
    # targets of targets of source node
    dataframe[dataframe[,source_column_number] %in% dataframe[ dataframe[,source_column_number] %in% source, target_column_number], 
              target_column_number]
  )
}
# function_three_intermediate_nodes_central_node("TP53","ATM",pypath_dataframe)
# this outputs 2 central nodes of 2 paths

# find the path with 3 intermediate nodes, by providing the source and the target
# the function outputs the row index of interactions leading from the source to the target
function_three_intermediate_nodes_path <- function(source, target, dataframe, source_column_number, target_column_number) {
  
  if (sum(class(dataframe)!="data.frame")>0){
    dataframe <- as.data.frame(dataframe)
  }
  
  if(missing(target_column_number)){
    target_column_number<-3
  }
  if(missing(source_column_number)){
    source_column_number<-1; 
  }
  
  central_node <- function_three_intermediate_nodes_central_node(source,target,dataframe)
  
  source_first_neighbor <- intersect( dataframe[dataframe[,source_column_number] %in% source, target_column_number],
                                      dataframe[dataframe[,target_column_number] %in% central_node, source_column_number])
  target_first_neighbor <- intersect( dataframe[dataframe[,target_column_number] %in% target, source_column_number], 
                                      dataframe[dataframe[,source_column_number] %in% central_node, target_column_number])
  # source -> source's 1st neighbour row
  c( # source -> source's 1st neighbour row
    which( dataframe[,source_column_number] %in% source & dataframe[, target_column_number] %in% source_first_neighbor ),
    # source's 1st neighbour row -> central node
    which( dataframe[,source_column_number] %in% source_first_neighbor & dataframe[,target_column_number] %in% central_node ),
    # central node -> target's first neigbour
    which(dataframe[,source_column_number] %in% central_node & dataframe[,target_column_number] %in% target_first_neighbor),
    # target's first neigbour -> target
    which( dataframe[,source_column_number] %in% target_first_neighbor & dataframe[,target_column_number] %in% target ) )
}


#####################
#####################

# input is the output of the function function_two_intermediate_nodes_path 
function_threestep_paths_table <- function(threestep_paths_sif,source,target,colnames_table){
  # pair rows with source nodes with their targets
  x <- threestep_paths_sif[threestep_paths_sif$source_hgnc %in% source,]; colnames(x)<-paste0(colnames(threestep_paths_sif),1,sep="")
  y <- threestep_paths_sif; colnames(y) <- paste0(colnames(threestep_paths_sif),2,sep="")
  z <- merge(x,y,by.x="target_hgnc1",by.y="source_hgnc2",sort=F)[,c(colnames(x),colnames(y)[!grepl("source",colnames(y))])]
  # pair this 3 element table to 
  df_target_rows <- threestep_paths_sif[threestep_paths_sif$target_hgnc %in% target,]
  colnames(df_target_rows) <- paste0(colnames(threestep_paths_sif),3,sep="")
  screen_model_threestep_paths <- merge(z,df_target_rows,by.x="target_hgnc2",by.y="source_hgnc3",sort=F)[,
    c(colnames(z),colnames(df_target_rows)[!grepl("source",colnames(df_target_rows))])]
  # order
  screen_model_threestep_paths <- screen_model_threestep_paths[order(screen_model_threestep_paths[,7],screen_model_threestep_paths[,1],
                                                                     screen_model_threestep_paths[,3], screen_model_threestep_paths[,5]),]
  # remove loops, like x1 -> x2 -> x1 -> x3
  
  if (prod(dim(screen_model_threestep_paths))>0){
  screen_model_threestep_paths <- screen_model_threestep_paths[!apply(apply(screen_model_threestep_paths[,
                    seq(1,ncol(screen_model_threestep_paths),2)][2:length(seq(1,ncol(screen_model_threestep_paths),2))],1,duplicated),2,sum),]
  } else {
    screen_model_threestep_paths <- as.data.frame(t(rep(NA,7)))
  }
  # rownames(screen_model_threestep_paths) <- c()
  colnames(screen_model_threestep_paths) <- rev(rev(colnames_table)[1:length(colnames(screen_model_threestep_paths))])
  screen_model_threestep_paths
}

######################
#####################

# function to generate all full paths involving three intermediate nodes, ie. of 4 steps

function_fourstep_paths_table <- function(screen_model_fourstep_paths_sif, source, target, colnames_table){
# x0 -> x1
df_step1 <- screen_model_fourstep_paths_sif[screen_model_fourstep_paths_sif$source_hgnc %in% source,]
colnames(df_step1) <- paste0(colnames(df_step1),1,sep="")
# x3 -> M(odel gene)
df_step4 <- screen_model_fourstep_paths_sif[screen_model_fourstep_paths_sif$target_hgnc %in% target,]
colnames(df_step4) <- paste0(colnames(screen_model_fourstep_paths_sif),4,sep="")
# x1 -> x2
df_step2 <- screen_model_fourstep_paths_sif[screen_model_fourstep_paths_sif$source_hgnc %in% df_step1$target_hgnc1,]
colnames(df_step2) <- paste0(colnames(screen_model_fourstep_paths_sif),2,sep="")
# x2 -> x3
df_step3 <- screen_model_fourstep_paths_sif[screen_model_fourstep_paths_sif$source_hgnc %in% df_step2$target_hgnc2 &
                                              screen_model_fourstep_paths_sif$target_hgnc %in% df_step4$source_hgnc4,]
colnames(df_step3) <- paste0(colnames(screen_model_fourstep_paths_sif),3,sep="")

df_step1_step2 <- merge(df_step1,df_step2,
                    by.x=colnames(df_step1)[target_column_number],
                    by.y=colnames(df_step2)[source_column_number],sort=F)[,c(colnames(df_step1),colnames(df_step2)[!grepl("source",colnames(df_step2))])]

df_step3_step4 <- merge(df_step3,df_step4,
                    by.x=colnames(df_step3)[target_column_number],
                    by.y=colnames(df_step4)[source_column_number],sort=F)[,c(colnames(df_step3),colnames(df_step4)[!grepl("source",colnames(df_step4))])]

screen_model_fourstep_paths <- merge(df_step1_step2,df_step3_step4,by.x=colnames(df_step1_step2)[length(df_step1_step2)],
      by.y=colnames(df_step3_step4)[source_column_number],sort=F)[,c(colnames(df_step1_step2),
                                                                     colnames(df_step3_step4)[!grepl("source",colnames(df_step3_step4))])]

if (prod(dim(screen_model_fourstep_paths))>0){
  # remove loops (like x4->x3->x4->x2->M)
  screen_model_fourstep_paths <- screen_model_fourstep_paths[!apply(
    apply(screen_model_fourstep_paths[,seq(1,ncol(screen_model_fourstep_paths),2)[2:length(seq(1,ncol(screen_model_fourstep_paths),2))]],
          1,duplicated),2,sum),]
  } else {
    screen_model_fourstep_paths <- as.data.frame(t(rep(NA,9)))
}

screen_model_fourstep_paths <- screen_model_fourstep_paths[order(screen_model_fourstep_paths[,ncol(screen_model_fourstep_paths)],screen_model_fourstep_paths[,1], 
                                  screen_model_fourstep_paths[,3],screen_model_fourstep_paths[,5],screen_model_fourstep_paths[,7]),]
colnames(screen_model_fourstep_paths) <- rev(rev(colnames_table)[1:length(colnames(screen_model_fourstep_paths))])
rownames(screen_model_fourstep_paths) <- c()
screen_model_fourstep_paths
}

#########################

# create concatenated table for screen-screen connections

function_internal_paths_model_nodes <- function(full_network,source, target,source_column_number,target_column_number, colnames_table){

# one intermediate (two step) paths
model_model_twostep_paths <- function_table_oneintermediate_paths(source, target,full_network,source_column_number,target_column_number, colnames_table)
# two intermed (three step) paths
model_model_threestep_paths_sif <- full_network[function_two_intermediate_nodes_path(source, target, full_network,source_column_number,target_column_number),
                                              source_column_number:target_column_number]
model_model_threestep_paths_sif <- model_model_threestep_paths_sif[!duplicated(model_model_threestep_paths_sif),]
# turn this into a table: M -> x1 -> x2 -> M
model_model_threestep_paths <- function_threestep_paths_table(model_model_threestep_paths_sif, source, target, colnames_table)

# three intermed (four step) paths
model_model_fourstep_paths_sif <- full_network[function_three_intermediate_nodes_path(source, target, full_network,source_column_number,target_column_number),
                                             source_column_number:target_column_number]
model_model_fourstep_paths_sif <- model_model_fourstep_paths_sif[!duplicated(model_model_fourstep_paths_sif),]
# turn this into a table: M -> x1 -> x2 -> M
model_model_fourstep_paths <- function_fourstep_paths_table(model_model_fourstep_paths_sif, source, target, colnames_table)

# bind dataframes of internal paths (bind_rows is dplyr package)
model_model_all_paths <- bind_rows(model_model_twostep_paths, model_model_threestep_paths,model_model_fourstep_paths)[,colnames(model_model_fourstep_paths)]
model_model_all_paths

}

######################

function_screen_model_paths <- function(full_network,source,target,source_column_number,target_column_number, colnames_table){

# Direct connections (one step paths, first neighbors)
screen_model_direct_connections <- full_network[as.data.frame(full_network)[,source_column_number] %in% source &
                                      as.data.frame(full_network)[,target_column_number] %in% target, source_column_number:target_column_number]
screen_model_direct_connections <- screen_model_direct_connections[!duplicated(screen_model_direct_connections),]
screen_model_direct_connections <- screen_model_direct_connections[order(function_collapse_df_column_to_vector(screen_model_direct_connections[,source_column_number]), 
                                             function_collapse_df_column_to_vector(screen_model_direct_connections[,target_column_number])),]
# uniform column names
colnames(screen_model_direct_connections) <- rev(rev(colnames_table)[source_column_number:length(colnames(screen_model_direct_connections))])

#########################
#########################
# Connections with one intermediate node (two step paths/2nd neighbors) connecting screen genes to model genes
screen_model_twostep_paths <- function_table_oneintermediate_paths(source, target,
                                                                   full_network,source_column_number,target_column_number, colnames_table)

#########################
#########################
# Connections with two intermediate nodes (three step paths/3rd neighbors) connecting screen genes to model genes

# this outputs the rows that constitute the paths, eg.: x0(screen gene)->x1, x1->x2; x2->x3 (model gene)
screen_model_threestep_paths_sif <- full_network[function_two_intermediate_nodes_path(source, target, full_network,source_column_number,target_column_number),
                                               source_column_number:target_column_number]
screen_model_threestep_paths_sif <- screen_model_threestep_paths_sif[!duplicated(screen_model_threestep_paths_sif),]
# turn this into a table: x0 -> x1 -> x2 -> x3
screen_model_threestep_paths <- function_threestep_paths_table(screen_model_threestep_paths_sif, source, target, colnames_table)

#########################
#########################
# Connections with three intermediate node (four step paths/4th neighbors) connecting screen genes to model genes

# source <- source; target <- target
screen_model_fourstep_paths_sif <- full_network[function_three_intermediate_nodes_path(source,target,full_network,source_column_number,target_column_number),
                                              source_column_number:target_column_number]
screen_model_fourstep_paths_sif <- screen_model_fourstep_paths_sif[!duplicated(screen_model_fourstep_paths_sif),]
# turn this into a table: x0 -> x1 -> x2 -> M
screen_model_fourstep_paths <- function_fourstep_paths_table(screen_model_fourstep_paths_sif, source, target, colnames_table)

#########################
# rbind all screen-model tables
screen_model_all_paths <- bind_rows(screen_model_direct_connections, screen_model_twostep_paths,
                                    screen_model_threestep_paths,screen_model_fourstep_paths)[,colnames(screen_model_fourstep_paths)]
screen_model_all_paths
}
