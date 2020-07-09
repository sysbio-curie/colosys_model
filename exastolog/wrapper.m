%% ExaStoLog
% This file contains the commands to run the functions 
% calculating the stationary solution of stochastic logical models, 
% plot results and perform sensitivity analysis (wrt transition rates)
%
% READ the tutorial at: https://github.com/mbkoltai/exact-stoch-log-mod
% author: Mihaly Koltai, mihaly.koltai@curie.fr

% go to the folder of the file
editor_service=com.mathworks.mlservices.MLEditorServices; editor_app=editor_service.getEditorApplication;
active_editor=editor_app.getActiveEditor; storage_location=active_editor.getStorageLocation;
file=char(storage_location.getFile); path_to_toolbox=fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS to PATH: required toolboxes are in <toolboxes.zip>
add_toolboxes_paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model set-up

% models can be defined 
% A) by entering the list of nodes and their
% corresponding rules as a cell of strings, using MATLAB logical notation ('&', '|', '~', '(', ')'),
% for instance a toy model with cyclic attractor:
% nodes={'A','B','C'}; % rules={'~A','A','A&~C'}
% nodes={'A','B','C','D'}; rules={'~B','~A&C','B|C&~D','C'}

% LIST of MODELS
model_name_list={'mammalian_cc', ...
'breast_cancer_zanudo2017',....
'EMT_cohen_ModNet',...
'sahin_breast_cancer_refined',...
'krasmodel15vars'}; % 
% name of the model
model_index=numel(model_name_list);
model_name=model_name_list{model_index};
% 'dnarepair_rodriguez_15nodes',...

% where to save figures
plot_save_folder=strcat('doc/sample_plots/',model_name,'/');

% model read in from an existing BOOLNET file
[nodes,rules]=fcn_bnet_readin(strcat('model_files/',model_name,'.bnet')); 

% once we have list of nodes and their logical rules, check if all variables referred to by rules found in list of nodes:
fcn_nodes_rules_cmp(nodes,rules)

% if yes, we generate a function file, which will create the model
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)

% from the model we generate the STG table, that is independent of values of transition rates
% this takes ~5 seconds for 20 node model (but need to do it only once)
tic; stg_cell=fcn_build_stg_cell(truth_table_filename,nodes); toc

% number of transitions
log10(sum(sum(cellfun(@(x) numel(x),stg_cell))))
% density of transition matrix
sum(sum(cellfun(@(x) numel(x),stg_cell)))/(2^(2*numel(nodes)))
% visualize transition matrix
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)
% save: export_fig(strcat(save_folder,model_name,'_A_sparse.pdf'),'-transparent','-nocrop','-r350')

%% generate transition matrix from existing STG table

% to define transition rates, we can select given rates to have different values than 1, or from randomly chosen
% name of rates: 'u_nodename' or 'd_nodename'
% chosen_rates={'u_ERBB1','u_ERBB2','u_ERBB3'}; chosen_rates_vals=zeros(size(chosen_rates));
% OR leave them empty: 
chosen_rates=[]; chosen_rates_vals=[];

% then we generate the table of transition rates: first row is the 'up'rates, second row 'down' rates, in the order of 'nodes'
% ARGUMENTS
distr_type={'uniform','random'}; % <uniform> assigns a value of 1 to all params. other option: <random>
meanval=[]; sd_val=[]; % if 'random' is chosen, the mean and standard dev of a normal distrib has to be defined 
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);
% meanval=1; sd_val=1; transition_rates_table=fcn_trans_rates_table(nodes,'random',meanval,sd_val,chosen_rates,chosen_rates_vals)

% build transition matrix A with parameter values
tic; [A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,''); toc
% if we want the kinetic matrix too, this is the 2nd output of the function
% tic; [A_sparse,K_sparse]=fcn_build_trans_matr(stg_table,transition_rates_table,'kinetic'); toc
    
% VISUALIZE transition matrix
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)
% fcn_save_fig('A_sparse',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

size_limit_mb=1; fcn_objects_memory_size(whos,size_limit_mb)

%% define initial condition

n_nodes=numel(nodes);
% truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);

% define some nodes with a fixed value and a probability <dom_prob>: 
% states satisfying this condition will have a total initial probability of <dom_prob>, the other states 1-dom_prob

% CELL CYCLE MODEL
% initial state specifying all variables:
% initial_fixed_nodes={'CycD','Rb_b1','Rb_b2','p27_b1','p27_b2','Cdh1','Skp2','E2F','CycE','CycA','CycB','Cdc20','UbcH10'}; 
% initial_fixed_nodes_vals=[ones(1,7) zeros(1,6)];
% 
% initial state specifying some variables:: CycE=0 & CycA=0 & CycB=0 & Cdh1=1 & % Rb=1 & p27=1
% initial_fixed_nodes={'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}; initial_fixed_nodes_vals=[0 0 0 1 1 1 1 1];
%
%
% KRAS model, WT: {'cc','KRAS'}, [1 0]. Mutant: {'cc','KRAS'}, [1 0]
% KRAS 15 nodes mutant
% initial_fixed_nodes={'cc','KRAS','cell_death'}; initial_fixed_nodes_vals=[1 1 0];

initial_fixed_nodes_list={ {'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}, ... % mammalian_cc
 {'Alpelisib', 'Everolimus','PIM','Proliferation','Apoptosis'},...  % breast_cancer_zanudo2017
 {'ECMicroenv','DNAdamage','Metastasis','Migration','Invasion','EMT','Apoptosis','Notch_pthw','p53'}, ... % EMT_cohen_ModNet 
 {'EGF','ERBB1','ERBB2','ERBB3','p21','p27'}}; % sahin_breast_cancer_refined % 'EGF','ERBB1','ERBB2','ERBB3','ERBB_12','ERBB_13','ERBB_23','CDK6','p21','p27'

initial_fixed_nodes_vals_list={[0 0 0 1 1 1 1 1], ...   % mammalian_cc
     [0 1 0 zeros(1,2)],...  % breast_cancer_zanudo2017
    [1 1 zeros(1,5) 1 0],... % EMT-Cohen model: [0/1 0/1 zeros(1,5)]
    [1 0 0 0 1 1]}; % 1 zeros(1,numel(initial_fixed_nodes_list{model_index})-3) 1 1
initial_fixed_nodes=initial_fixed_nodes_list{model_index}; initial_fixed_nodes_vals=initial_fixed_nodes_vals_list{model_index};

% what is the probability of this state, (eg. dom_prob=0.8, ie. 80% probability)
dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'uniform','random'};
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to the states with the fixed nodes having the defined values
tic; x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_types{1},plot_flag); toc

% completely random initial condition:
% x0=zeros(2^n_nodes,1); x0=rand(2^numel(nodes),1); x0=x0/sum(x0);
% completely uniform initial condition
% x0=ones(2^numel(nodes),1)/(2^numel(nodes));

%% CALCULATE STATIONARY STATE

% get the subnetworks of the STG and topologically sort them
% ARGUMENTS:
% transition matrix: A
% table of transition rates: transition_rates_table
% initial conditions: x0
tic; stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0); toc
% <stg_sorting_cell> contains
% {subnetws: which subgraph each state belongs to, #1
% scc_submat_cell: states in the subgraphs, #2
% nonempty_subgraphs: which subgraphs are populated by the initial condition, #3
% sorted_vertices_cell: states (vertices) topologically sorted in each nonempty subgraph, #4
% cyclic_sorted_subgraphs_cell: sorted states within cycles} #5

%% calculate stationary solution
tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0); toc

% OUTPUTS
% stat_sol: stationary solution for all the states
% term_verts_cell: index of nonzero states. If STG is disconnected the nonzero states in these disconn subgraphs are in separate cells
% cell_subgraphs: indices of states belonging to disconnected subgraphs (if any)

% query size of objects larger than x (in Mbytes)
% objects_mem=whos; size_limit_mb=1; fcn_objects_memory_size(objects_mem,size_limit_mb)

% sum the probabilities of nonzero states by nodes, both for the initial condition and the stationary solution
% ARGUMENTS
% initial conditions: x0
% stat_sol: stationary solution for all the states
% nodes: list of nodes
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0');

% Checked with MaBoSS simuls, results are identical (up to 1% deviation).
% Comparing with simulation of mammalian cell cycle model with 12 nodes: 
% look in folder 'doc/sample_plots/maboss'

% for dynamical simulation of small (<5) models use:
% [~,x]=ode45(@(t,x)fcn_K_ode(t,x,full(K_sparse)),tspan,x0);

%% PLOTTING RESULTS

% this function plots 2 or 3 subplots:
% 1) transition matrix (optional)
% 2) stationary probability of the model's states
% 3) stationary probability of the model's variables (having the value of 1)

% PLOT A/K and stat solutions
% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% prob_thresh: minimal value for probability to display (if non-empty, only plot nonzero states, useful for visibility if many states)
sel_nodes=[];
min_max_col=[0 1]; barwidth_states_val=0.8;fontsize=[24 40 20]; % [fontsize of plot, fontsize of titles, fontsize of binary states]
plot_settings = [fontsize barwidth_states_val min_max_col]; prob_thresh=0.02;
% WARNING!!! if more than 12 nodes, generating the figure for A/K can be time-consuming
% leave first variable (A_sparse) empty ([]) to have plot without matrix
figure('name','A_K_stat_sol')
fcn_plot_A_K_stat_sol(A_sparse,nodes,sel_nodes,stat_sol,x0,plot_settings,prob_thresh)

% SAVE
% enter any string for the last argument to overwrite existing plot!!
if exist(plot_save_folder,'dir')==0; mkdir(plot_save_folder); end
fig_file_type={'.png','.eps','.pdf','.jpg','.tif'}; 
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually
resolution_dpi='-r350'; % magnification=0.8;  strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
fcn_save_fig('single_solution_states_nodes_stat_sol',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

%% PLOT binary heatmap of nonzero stationary states with their probability

% ARGUMENTS of function:
% stat_sol: vector of stationary solutions
% prob_thresh: probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.01;  % []; % 0.05;
% term_verts_cell: index of subgraphs for stable states
% nodes: name of nodes
% sel_nodes: nodes to show. if none selected, then all nodes shown
sel_nodes=[]; % setdiff(2:numel(nodes)-1,[find(strcmp(nodes,{'Rb_b2'})) find(strcmp(nodes,{'p27_b2'}))]);
% plot_param_settings
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=26; fontsize=36; hor_gap=0.02; bottom_marg=0.31; left_marg=0.22; 
plot_param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% tight_subplot_flag: want to use tight subplot? | ranking_flag: order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';

% PLOT
figure('name','statsol_binary_heatmap')
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,...
                        term_verts_cell, nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag);

% SAVE
% if <overwrite_flag> non-empty then existing file with same name is overwritten. 
overwrite_flag='yes'; 
% for <resolution> you can enter dpi value manually
% magnification=0.8; resolution_dpi=strcat('-r',num2str(magnification*get(0, 'ScreenPixelsPerInch')));
resolution_dpi='-r350'; fcn_save_fig('binary_heatmap_states',plot_save_folder,fig_file_type{3},overwrite_flag,resolution_dpi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sensitivity analysis: one-dimensional parameter scans

% select the nodes whose parameters we want to scan in:
% all nodes that have actual transitions
% par_inds_table=unique(stg_table(:,[3 4]),'rows');
popul_subgraphs=cellfun(@(x) sum(ismember(find(x0>0),x)), cell_subgraphs)>0;
subgraph_states=cell2mat(cell_subgraphs(popul_subgraphs)');
% subgraph_states=cell2mat(cell_subgraphs(~cellfun(@(x) isempty(x),term_verts_cell))');
state_trans_freq=cell2mat(cellfun(@(x) sum(ismember(x,subgraph_states)), stg_cell','un',0));
[a,b,~]=find(state_trans_freq>0);
par_inds_table=[a,b];
% unique(stg_table(ismember(stg_table(:,1), subgraph_states) | ismember(stg_table(:,2), subgraph_states),3:4),'rows');

% most common transitions
% for k=1:size(par_inds_table,1)
%     param_freq(k) = sum(stg_table(:,3)==par_inds_table(k,1) & stg_table(:,4)==par_inds_table(k,2));
% end
% top n most frequent transitions
[~,top_freq_trans_rates]=sort(cell2mat(arrayfun(@(x) state_trans_freq(par_inds_table(x,1),par_inds_table(x,2)),...
                                1:size(par_inds_table,1),'un',0)),'descend');

% all rates that have corresponding transitions
scan_params=unique(par_inds_table(:,1))';
% all nodes that have transitions: unique(par_inds_table(:,1))'; 
% most frequent: scan_params=par_inds_table(top_freq_trans_rates(1:6),1)';
% Zanudo: find(ismember(nodes,{'AKT','SGK1','TSC','FOXO3','BIM','BAD','mTORC1','PI3K','PRAS40'})); 
% all nodes that have transitions, except phenotypes: setdiff(unique(par_inds_table(:,1))',find(ismember(nodes,{'Apoptosis','Proliferation'})))
% selected nodes: find(ismember(nodes,{'AKT','SGK1','TSC','FOXO3','BIM','BAD','mTORC1'})); % 
scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0); 
% how many transition rates we'll scan? sum(cellfun(@(x) numel(x),scan_params_up_down))
%
% num2cell(repelem([1 2],numel(scan_params),1),2)'; % both up and down rates
% num2cell(ones(1,numel(scan_params))); % only up 
% num2cell(2*ones(1,numel(scan_params))); % only down
% {[1 2], 1, [1 2]}; % manually selected

% top frequency trans rates:
% scan_params=par_inds_table(top_freq_trans_rates,1)'; 
% scan_params_up_down=arrayfun(@(x) par_inds_table(top_freq_trans_rates(par_inds_table(top_freq_trans_rates,1)==x),2)',scan_params,'un',0);

% min and max of range of values; resolution of the scan; linear or logarithmic sampling
parscan_min_max = [1e-2 1e2]; n_steps=5; sampling_types={'log','linear'}; 

% FUNCTION for generating matrix of ordered values for the parameters to scan in
% [scan_par_table,scan_par_inds,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,transition_rates_table);
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,nodes,sampling_types{1},parscan_min_max,n_steps);
% to set values in one/more columns manually (order is {3,1}->{3,2}->{4,1}->{4,2} etc)

% FUNCTION for 1-DIMENSIONAL PARSCAN
% ARGUMENTS
% stg table: generated above, state transition graph
% transition_rates_table:default values of transition rates
% initial conditions: x0
% stg_table: generated by stg_table=fcn_build_stg_table(truth_table_filename,nodes);
% [~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
tic;
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_cell,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);
toc;

%% PLOT RESULTS of 1-dimensional parameter scan on heatmap/lineplot BY PARAMETERS

%%% FIRST PLOT TYPE: each subplot is a lineplot of node or state values as a function of a parameter's value, 
%%% with a defined minimal variation

% index of nonzero states
nonzero_states_inds=find(stat_sol>0);
% plot parameters
% [0.06 0.03],[0.03 0.03],[0.02 0.01]
height_width_gap=[0.08 0.03]; bott_top_marg =[0.05 0.05]; left_right_marg=[0.04 0.01];
params_tight_subplots={height_width_gap bott_top_marg left_right_marg};
% [fontsize_axes,fontsize_title,legend_fontsize,linewidth,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={24,34,24,4,params_tight_subplots,model_name};
state_or_node_flags={'nodes','states'}; 
% cutoff for minimal variation to show a variable
diff_cutoff=0.15;
figure('name','onedim parscan by param')
[fig_filename,onedim_paramscan_output_cell]=fcn_onedim_parscan_plot_by_params(state_or_node_flags{1},...
                                      stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                      nonzero_states_inds,parscan_matrix,nodes,...
                                      scan_params,scan_params_up_down,... % selected parameters
                                      diff_cutoff,... % minimal variation for variable to be shown on plot
                                      plot_param_settings);
% SAVE figure
% resolution_dpi='-r350'; fcn_save_fig(strcat(fig_filename,'_r350'),plot_save_folder,fig_file_type{2},'overwrite',resolution_dpi);

%% PLOT RESULTS of 1-by-1 parameter scan on heatmap/lineplot BY VARIABLES

%%% SECOND PLOT TYPE: show stationary value/response coefficient of 1 variable (state or node) on 1 subplot, as a fcn of all relevant parameters
% nonzero states of the model
% nonzero_states=unique(cell2mat(stationary_state_inds_scan(:)'))';
nonzero_states_inds=find(stat_sol>0);
sensit_cutoff=0.1; % minimal value for response coefficient (local sensitivity) or for the variation of node/state values
% select parameters of plot
height_width_gap=[0.1 0.04]; bott_top_marg=[0.03 0.1]; left_right_marg=[0.07 0.02];
% [fontsize_axes,fontsize_title,params_tight_subplots(leave empty if not installed),model_name]
plot_param_settings={30,30,{height_width_gap bott_top_marg left_right_marg},model_name,'colorbar'}; 
% plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
% if want to loop through all plot types: all_opts_perm=[[1 1 1]; unique([perms([1 1 2]); perms([2 2 1])],'rows'); [2 2 2]];
plot_type_options=[1 1 1];
figure('name',strjoin(arrayfun(@(x) plot_types{x}{plot_type_options(x)}, 1:numel(plot_type_options), 'un',0),'_'));
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                                   stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                                   nonzero_states_inds,parscan_matrix,nodes,...
                                                   scan_params,scan_params_up_down,... % transition rates to scan in
                                                   sensit_cutoff,plot_param_settings);

% <resp_coeffs> dimensions: (parameters, values,nodes), so eg. resp_coeffs(:,:,7)=resp. coeff values across param ranges for the 7th node

% SAVE figure
% resolution_dpi='-r350'; 
% fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(sensit_cutoff),'.','p')),plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% Multidimensional param sampling at UNIFORM distances (2-dimensions in example below)

% # of cols in <all_par_vals_lhs> has to be same as # of elements in <scan_params_up_down>
% example: 2-dimensional uniform scan in u_Notch_pthw and u_p53
n_scanvals=5; scanvals=logspace(-2,2,n_scanvals);  % with zero [0 logspace(-2,2,n_scanvals-1)]
meshgrid_scanvals=meshgrid(scanvals,scanvals);

paramsample_table=[repelem(scanvals,n_scanvals)' reshape(reshape(repelem(scanvals,n_scanvals),n_scanvals,n_scanvals)',n_scanvals^2,1)]; 
multiscan_pars=[11 13]; multiscan_pars_up_down={1 1};

disp_var=5; % show at every n% the progress
[stat_sol_paramsample_table,stat_sol_states_paramsample_table]=fcn_calc_paramsample_table(paramsample_table,multiscan_pars,...
                                                                multiscan_pars_up_down,transition_rates_table,stg_cell,x0,disp_var);

% PLOT heatmap of selected variables 
sel_nodes=4; % [4 8 13 15];
% plot_settings: [fontsize on plot&axes, fontsize on axes, fontsize of subplot titles, axes tick fontsize]
plot_settings=[28 30 40]; figure('name','2D scan')
fcn_plot_twodim_parscan(stat_sol_paramsample_table,scanvals,multiscan_pars,multiscan_pars_up_down,nodes,sel_nodes,plot_settings)

% SAVE
% twodim_parscan_Apoptosis_Metastasis
resolution_dpi='-r350'; file_name_prefix=strcat('twodim_parscan_',strjoin(nodes(sel_nodes),'_'));
fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);

%% multidimensional parameter scan: LATIN HYPERCUBE SAMPLING (random multidimensional sampling within given parameter ranges)

% WHICH PARAMETERS to scan simultaneously?
% all relevant parameters: 
% scan_params=unique(stg_table(:,3))'; scan_params_up_down=num2cell(repelem([1 2],numel(scan_params),1),2)'; 
% scan_params_up_down_sensit=arrayfun(@(x) scan_params_up_down_sensit{x}', 1:numel(scan_params_up_down_sensit),'un',0)

% PERFORM Latin Hypercube sampling (LHS) SAMPLING
sampling_types={'lognorm','linear','logunif'}; sampling_type=sampling_types{3};
% par_min_mean: minimum or in case of lognormal the mean of distribution. Can be a scalar or a vector, 
% if we want different values for different parameters
% max_stdev: maximum or in case of lognormal the mean of distribution. 
% Can be a scalar or a vector, if we want different values for different parameters
%
% for 'lognorm' and 'logunif' provide the LOG10 value of desired mean/min and stdev/max!!, ie. -2 means a mean of 0.01
par_min_mean=-2; % repmat(1.5,1,numel(cell2mat(scan_params_up_down(:)))); par_min_mean(4)=3; 
max_stdev=2; % repmat(0.5,1,numel(cell2mat(scan_params_up_down(:))));
% <lhs_scan_dim>: number of param sets
lhs_scan_dim=1000;
[all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan]=... % outputs
    fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
                                            scan_params_sensit,scan_params_up_down_sensit, ... % transition rates
                                            transition_rates_table,stg_cell,x0,nodes);
                    
%% SCATTERPLOTS of STATE or NODE values as a function of the selected parameters, with the trendline shown (average value per parameter bin)

% which variable to plot?
var_ind=3;
% STATES or NODES? <scan_values>: values to be plotted
scan_values=stat_sol_states_lhs_parscan; % model variables: stat_sol_nodes_lhs_parscan; model states: stat_sol_states_lhs_parscan

% PLOT
sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
% file_name_prefix=strcat('LHS_parscan_scatterplot_trend_',nodes{var_ind}); 
file_name_prefix=strcat('LHS_parscan_scatterplot_trend_state',num2str(var_ind));
% param_settings: [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
param_settings = [50 6 24 size(stat_sol_states_lhs_parscan)];

figure('name',num2str(var_ind))
fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,...
        scan_params_sensit,scan_params_up_down_sensit,nodes,sampling_type,param_settings)

resolution_dpi='-r350'; fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);
% fcn_save_fig('LHS_parscan_trend_AKT2',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% calculating & plotting (heatmap) correlations between variables 

% ARGUMENTS
% sel_nodes: name of selected nodes (pls provide in ascending order) (if left empty, all shown)
sel_nodes=[3 7 8 10 11 13:15 17:20]; 
% plot_settings: [fontsize on plot, fontsize on axes/labels]
plot_settings=[NaN 26 32]; 
plot_type_flag={'var_var','heatmap'}; % this is plotting the heatmap of correlations between variables
figure('name',strjoin(plot_type_flag))
[varvar_corr_matr,p_matrix_vars]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,stat_sol_nodes_lhs_parscan,...
                                            nodes,sel_nodes,[],[],[],plot_settings);

fig_prefix=strcat(strjoin(plot_type_flag,'_'),'_corrs');
% resolution_dpi='-r350'; fcn_save_fig(fig_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi);


%% linear or lin-log regression of VARIABLES as fcn of PARAMETERS: VARIABLE=f(PARAMETER), the function plots R squared

plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
sel_nodes=[]; % setdiff(3:numel(nodes),[4:6 10:12 16]); %setdiff(3:numel(nodes),[3 5 6 7 9 11 12 16]); % 3:numel(nodes); % scan_params_sensit;
% plot_settings=[fontsize,maximum value for heatmap colors], if plot_settings(3)=NaN, then max color automatically selected
plot_settings=[30 30 0.29]; 
% if regression type is 'linlog', then the fit is y = a + b*log10(x)
regr_types={'log','linear'}; % log recommended if parameter values log-uniformly distributed in sampling
figure('name',strjoin(plot_type_flag))
scan_values=stat_sol_states_lhs_parscan; % stat_sol_nodes_lhs_parscan
[r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,all_par_vals_lhs,scan_values,...
                                 nodes,sel_nodes,... % which nodes
                                 scan_params_sensit,scan_params_up_down_sensit, ... % parameters (CAREFUL that they are same as in LHS!)
                                 regr_types{1},plot_settings);

% savefig(strcat(plot_save_folder,fig_prefix))                           
% fig_prefix=strjoin(plot_type_flag,'_'); resolution_dpi='-r350'; 
% fcn_save_fig(fig_prefix,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
   
%% Sobol total sensitivity metric                                            

% On Sobol total sensitivity index see: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
% This metric indicates how much of the total variance in a variable is due to variation in a given parameter
% We calculate here the usual numerical approximation of analytical equivalent from Monte Carlo sampling.
% From the LHS sampling above we take the matrices of parameter sets and variable values:
% [parameter sets, variable values]: [all_par_vals_lhs,stat_sol_nodes_lhs_parscan]

% select only parameters with an R-squared over some value
% careful to use same parameters as for R^2 calculations!!!
[par_ind_table,sequential_indices_lhs,~] = fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);
% R^2>=0.1
r_sq_thresh=0.05;
par_ind_table_filtered=par_ind_table(sum(r_squared>r_sq_thresh)>0,:);
scan_params_filtered=unique(par_ind_table_filtered(:,1))'; % scan_params_unique=unique(par_ind_table(sum(r_squared>0.2)>0,1))'; 
scan_params_up_down_filtered=arrayfun(@(x) par_ind_table_filtered(par_ind_table_filtered(:,1)==x,2)', scan_params_filtered,'un',0);

% Sobol total sensitivity: calculated for one variable at a time
% selected nodes to display
sel_nodes=[]; % 3:numel(nodes); % scan_params_sensit;
sample_size=[]; % if left empty, the sample size is half of the original param scan <all_par_vals_lhs>
% how often should be progress displayed?
disp_freq=10;
plot_settings=[];
var_types={'node','state'}; % analysis for states or nodes
% to calculate Sobol total sensitivity we need <sample_size*numel(scan_params_up_down)> evaluations of the model
figure('name','sobol sensitivity index')
sobol_sensit_index=fcn_multidim_parscan_sobol_sensit_index([],var_types{2},...
                      all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
                      sample_size,... % # of calculations per parameter
                      sequential_indices_lhs,... % this is indices of transition rates in the original LHS
                      scan_params_filtered,scan_params_up_down_filtered,...% scan_params_sensit,scan_params_up_down_sensit
                      stg_cell,transition_rates_table,x0,nodes,sel_nodes,plot_settings,disp_freq);

% if already calculated <sobol_sensit_index> and only want to plot results, provide <sobol_sensit_index> as FIRST argument 
% PLOT SETTINGS: [fontsize_plot,fontsize_axes,fontsize_title, min_color(optional), max_color(opt), angle of x-axis labels];
plot_settings=[30 30 40 0 0.5 90];
fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_types{2},all_par_vals_lhs,[],[],[],...
                       sequential_indices_lhs,scan_params_filtered,scan_params_up_down_filtered,...
                       [],[],[],nodes,sel_nodes,plot_settings,[]);
xticklabels({'Metastasis','Apoptosis (p53)','Apoptosis (p63_73)'})

% SAVE
% resolution_dpi='-r350'; fcn_save_fig('sobol_sensitivity_index',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% PARAMETER FITTING: define data, initial guess

% scan_params_sensit=[11 13 15 16]; scan_params_up_down_sensit={2,1,2,[1 2]};
% names of selected transition rates
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% define data vector (generate some data OR load from elsewhere)
data_param_vals=lognrnd(1,1,1,numel(predictor_names)); 
% initial guess for parameters
init_par_vals=data_param_vals.*lognrnd(1,2,size(predictor_names)); 

% initial true value of variables/states, initial guess
var_type_flag='states'; % 'vars' 'states'
[y_data,y_init_pred,init_error]=fcn_param_fitting_data_initguess_error(var_type_flag,...
                                        x0,stg_cell,data_param_vals,init_par_vals,...
                                        stg_sorting_cell,nodes,predictor_names);

%% function handles for fitting
[fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_handles_fitting(var_type_flag,...
                                y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);

%% FITTING by simul anneal

% default values for fitting hyperparameters:
% struct('CoolSched',@(T) (0.8*T), 'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
% 'InitTemp',1,'MaxConsRej',1000, 'MaxSuccess',20,...
% 'MaxTries',300, 'StopTemp',1e-8, 'StopVal',-Inf, 'Verbosity',1);
fitting_arguments=struct('Verbosity',2, 'StopVal', init_error/10, 'MaxTries',30,'MaxConsRej',100);
% FIT
tic; [optim_par_vals,best_error,T_loss]=anneal(fcn_statsol_sum_sq_dev,init_par_vals,fitting_arguments); toc;

% RESULTS
% values with fitted parameters
[y_optim_param,~,~]=fcn_param_fitting_data_initguess_error(var_type_flag,x0,stg_cell,data_param_vals,optim_par_vals,...
                                            stg_sorting_cell,nodes,predictor_names);

% values of fitted variables: [initial guess, true values (data), fitted values]
% careful with dimensions: for states, these are column vectors
data_init_optim=[y_init_pred; y_data; y_optim_param]; 
% parameters: initial guess, true values, fitted values
param_sets=[init_par_vals;data_param_vals;optim_par_vals];

%% PLOT: simulated annealing

figure('name','param fitting (simul.ann.)'); 
% select nodes to plot (here we selected nodes that are not always 0 or 1)
sel_nodes=find(sum(data_init_optim)>0 & sum(data_init_optim)<3);
% PLOT fitting process
thres_ind=size(T_loss,1); % thres_ind=find(T_loss(:,2)<1e-2,1); 
plot_settings=[24 30]; 
% var_type_flag='vars'; % 'states'
figure('name','simul anneal')
fcn_plot_paramfitting(var_type_flag,data_init_optim,T_loss,nodes,sel_nodes,[1 2],thres_ind,plot_settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitting by initial numerical gradient

error_thresh_fraction=0.1; 	% what fraction of the initial error to stop?
step_thresh=[]; 	% what step # to stop? you can leave this empty 
% init_error_table: changes to initial error when increasing or decreasing parameter values
init_error_table=[]; % if we have it from previous fitting than feed it to fcn
% incr_resol_init: initial % change from initial param values to calculate the numerical gradient
% incr_resol: change in param values during gradient descent
incr_resol_init=0.15; incr_resol=0.03;

% FIT
% var_type_flag: 'states' or 'vars'
% careful that string and data type are consistent!!
[init_error_table,optim_pars_conv,statsol_parscan,error_conv]=fcn_num_grad_descent(var_type_flag,init_error_table,...
	{y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names},data_param_vals,...
	init_par_vals,incr_resol,incr_resol_init,error_thresh_fraction,step_thresh);

%% PLOT
% which vars/states to show, if empty all are shown
sel_nodes=[]; plot_settings=[24 30]; 
% if its states you fitted, take the transpose of ydata
data_init_optim=[statsol_parscan([1 end],:); y_data']; 
figure('name','numer grad_desc') % state_var_flags={'state','var'};
fcn_plot_paramfitting(var_type_flag,data_init_optim,error_conv,nodes,sel_nodes,[],[],plot_settings)

%% SAVE
fig_name=strcat('grad_descent_',var_type_flag,'_',num2str(numel(predictor_names)),'fittingpars');
fcn_save_fig(fig_name,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi)
