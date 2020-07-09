% commands for tutorial
% go to the folder of the file
editor_service=com.mathworks.mlservices.MLEditorServices; editor_app=editor_service.getEditorApplication;
active_editor=editor_app.getActiveEditor; storage_location=active_editor.getStorageLocation;
file=char(storage_location.getFile); path_to_toolbox=fileparts(file); cd(path_to_toolbox);

% FIRST: unzip the file 'toolboxes.zip' for external libraries
if exist('toolboxes.zip','file')>0; unzip toolboxes.zip; end

% go to the folder of the file
% editor_service=com.mathworks.mlservices.MLEditorServices; editor_app=editor_service.getEditorApplication;
% active_editor=editor_app.getActiveEditor; storage_location=active_editor.getStorageLocation;
% file=char(storage_location.getFile); path_to_toolbox=fileparts(file); cd(path_to_toolbox);

% ADD FUNCTIONS and external libraries (from 'toolboxes.zip') to PATH
add_toolboxes_paths

%% READ IN model

% names of models
model_name_list = {'mammalian_cc', ...
'breast_cancer_zanudo2017'....
'EMT_cohen_ModNet',...
'sahin_breast_cancer_refined',...
'krasmodel15vars',...
'krasreduced',...
'krasreduced_dnarep',...
'krasreduced_dnarep_pruned',...
'krasreduced_dnarep_simplif'}; %
% select the index of one model
model_index=9;
model_name=model_name_list{model_index};

% read in model from BOOLNET file
[nodes,rules]=fcn_bnet_readin(strcat('model_files/',model_name,'.bnet'));
% Compare if rules and nodes are consistent
fcn_nodes_rules_cmp(nodes,rules)

% where to save figures
plot_save_folder=strcat('doc/sample_plots/',model_name,'/');
% write file with logical rules
truth_table_filename='fcn_truthtable.m'; fcn_write_logicrules(nodes,rules,truth_table_filename)
% build STG
tic; stg_cell=fcn_build_stg_cell(truth_table_filename,nodes); toc
% density of STG: 
sum(sum(cellfun(@(x) numel(x),stg_cell)))/(2^(2*numel(nodes)))

%% choose transition rates

chosen_rates=[]; chosen_rates_vals=[];
% ARGUMENTS
% <uniform> assigns a value of 1 to all params. <random> samples from a lognormal distribution
distr_type={'uniform','random'};
% if <random> is chosen, the mean and standard dev of a normal distrib has to be defined
meanval=[]; sd_val=[];
transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);

% KRAS mutant?
kras_mut=0;
if kras_mut==1
    transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,{'d_KRAS','u_TP53','d_DSB_SSB','u_DNA_rep'},[0,1,1,1]);
    % transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,{'d_DSB_SSB','u_DNA_rep'},[1e1,1e1]);
end

%% BUILD transition matrix

tic; [A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,''); toc

% visualize
% spy(A_sparse); xlabel('model states'); ylabel('model states'); set(gca,'FontSize',24)

%% Defining initial conditions

% selected nodes for inital conditions
initial_fixed_nodes_list=...
{{'CycE','CycA','CycB','Cdh1','Rb_b1','Rb_b2','p27_b1','p27_b2'}, ... % mammalian_cc
{'Alpelisib', 'Everolimus','PIM','Proliferation','Apoptosis'},...  % breast_cancer_zanudo2017
{'ECMicroenv','DNAdamage','Metastasis','Migration','Invasion','EMT','Apoptosis','Notch_pthw','p53'}, ... % EMT
{'EGF','ERBB1','ERBB2','ERBB3','p21','p27'}, ... % sahin_breast_cancer_refined
{'cc','KRAS','DSB','CDC25B','g2m_trans','cell_death'},... % KRAS 15
{'KRAS','EGF','TGFb','GF','ATM_ATR','TP53','CDC25B_C','CDC25A','Proliferation','CASP3'},...
{'EGFR','RAS','DSB_SSB','TGFb','ATM_ATR','DNArepair','CDC25B_C','CDC25A','Proliferation','CASP3'},...
{'EGFR','RAS','DSB_SSB','TGFb','ATM_ATR','DNArepair','CDC25B_C','CDC25A','Proliferation','CASP3'},...
{'EGFR','TGFb','DSB_SSB','ATM_ATR','DNArepair','CDC25B_C','CDC25A','Proliferation','CASP3'}};

% values for selected nodes
initial_fixed_nodes_vals_list = {[0 0 0 1 1 1 1 1], ... % mammalian_cc
            [0 1 0 zeros(1,2)],...  % breast_cancer_zanudo2017
            [1 1 zeros(1,5) 1 0],... % EMT-Cohen model: [0/1 0/1 zeros(1,5)]
            [1 0 0 0 1 1],... % 1 zeros(1,numel(initial_fixed_nodes_list{model_index})-3) 1 1
            [1 1 1 0 0 0],...
            [1 1 1 1 0 0 0 0 0 0 0],...
            [1,1,1,1,0,0,0,0,0,0],...
            [1,1,1,1,0,0,0,0,0,0],...
            [1,1,1,0,0,0,0,0,0]};

% select the initial condition for the model we are working on
initial_fixed_nodes=initial_fixed_nodes_list{model_index}; 
initial_fixed_nodes_vals=initial_fixed_nodes_vals_list{model_index};

dom_prob=1;
% if <random> the probability is randomly distributed among states, if <uniform> uniformly
distrib_types={'random','uniform'}; n_distrib_type=2;
% if plot_flag non-empty, we get a bar plot of initial values
plot_flag='';
% function assigns a probability of <dom_prob> to states with the fixed nodes
x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,...
				dom_prob,nodes,distrib_types{n_distrib_type},plot_flag);

%% Topological sorting of STG and identification of cycles

stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);

% check size of created objects
size_limit_mb=1; % fcn_objects_memory_size(whos,size_limit_mb)

%% calculate stationary solution

tic; [stat_sol,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0); toc
% by model variables
[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0');

%% visualize solution

% ARGUMENTS
% matrix_input: [], K_sparse or A_sparse (kinetic or transition matrix)
% min_max_col: minimum and max color for heatmap
% fontsize: [font size on the heatmap, title font size for stationary solutions]
% barwidth_states_val: width of the bars for bar plot of stationary solutions of states
% sel_nodes: nodes to show. If left empty, all nodes are shown
% prob_thresh: minimal value for probability to display 
% (useful for visibility if many attractor states or large cyclic attractor(s))

% Call the function by:
sel_nodes=[];
min_max_col=[0 1]; barwidth_states_val=0.8;
% fontsize: [fontsize of plot, fontsize of titles, fontsize of binary states]
fontsize=[18 20 16];
plot_settings = [fontsize barwidth_states_val min_max_col]; prob_thresh=0.03;

figure('name','A_K_stat_sol')
fcn_plot_A_K_stat_sol(A_sparse,nodes,sel_nodes,stat_sol,x0,plot_settings,prob_thresh)

%%% SAVE figure
save_flag=0;
if save_flag==1
    if exist(plot_save_folder,'dir')==0; mkdir(plot_save_folder); end
    fig_file_type={'.png','.eps','.pdf','.jpg','.tif'};
% if <overwrite_flag> non-empty then existing file with same name is overwritten.
    overwrite_flag='yes';

% resolution of the figures (dpi)
    resolution_dpi='-r350';
% SAVE
    fcn_save_fig('single_solution_states_nodes_stat_sol',plot_save_folder,fig_file_type{1},overwrite_flag,resolution_dpi)
end

% to export to PDF you need to have GhostScript installed, install from: https://www.ghostscript.com/
% to export to eps requires pdftops, part of the Xpdf package, install from: % http://www.xpdfreader.com
% see also tutorial of export_fig at: https://github.com/altmany/export_fig/blob/master/README.md

%% Visualize binary heatmap of nonzero stationary states

% stat_sol: vector of stationary solutions
% prob_thresh: probability threshold for states to show (if left empty, all states shown)
prob_thresh=0.002;
% term_verts_cell: index of subgraphs for stable states
% nodes: name of nodes
% sel_nodes: nodes to show. if none selected, all nodes will be shown
sel_nodes=[];
% plot_param_settings
% num_size_plot: font size of 0/1s on the heatmap
% hor_gap: horizontal gap between terminal SCCs, bottom_marg: bottom margin, left_marg: left margin
numsize_plot=15; fontsize=22; hor_gap=0.02; bottom_marg=0.2; left_marg=0.06;
plot_param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
% tight_subplot_flag: want to use tight subplot? | ranking_flag: order states by probability?
tight_subplot_flag='yes'; ranking_flag='yes';

% PLOT
figure('name','statsol_binary_heatmap')
statsol_binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,...
term_verts_cell,nodes,sel_nodes,plot_param_settings,tight_subplot_flag,ranking_flag);

% SAVE
save_flag=0;
if save_flag==1
    disp('saving')
resolution_dpi='-r350';
fcn_save_fig('binary_heatmap_states',plot_save_folder,fig_file_type{1},overwrite_flag,resolution_dpi);
end

%% select transition rates for 1-dimensional parameter scan

popul_subgraphs=cellfun(@(x) sum(ismember(find(x0>0),x)), cell_subgraphs)>0;
subgraph_states=cell2mat(cell_subgraphs(popul_subgraphs)');
% subgraph_states=cell2mat(cell_subgraphs(~cellfun(@(x) isempty(x),term_verts_cell))');
state_trans_freq=cell2mat(cellfun(@(x) sum(ismember(x,subgraph_states)), stg_cell','un',0));
[a,b,~]=find(state_trans_freq>0); par_inds_table=[a,b];

% most common transitions
[~,top_freq_trans_rates]=sort(cell2mat(arrayfun(@(x) state_trans_freq(par_inds_table(x,1),par_inds_table(x,2)),...
                                1:size(par_inds_table,1),'un',0)),'descend');
% all
scan_params=unique(par_inds_table(:,1))';
% most common: scan_params=par_inds_table(top_freq_trans_rates(1:6),1)';
% by name: scan_params=find(ismember(nodes,{'Notch_pthw','p53','EMTreg','FOXO3','p63_73'}));

% up and down rates (all)
scan_params_up_down=arrayfun(@(x) par_inds_table(par_inds_table(:,1)==x,2)', scan_params,'un',0);

%% CALCULATE 1-dimensional scan

% min and max of range of values; resolution of the scan; linear or logarithmic sampling
parscan_min_max=[1e-2,1e2]; n_steps=3; sampling_types={'log','linear'};

% matrix of parameter values
parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,...
		nodes,sampling_types{1},parscan_min_max,n_steps);

% calculation
[stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,~]=... % stationary_state_inds_scan
    fcn_onedim_parscan_calc(stg_cell,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down);

%% PLOT 1-dimensional scan grouped by transition rates

% index of nonzero states
nonzero_states_inds=find(stat_sol>0);
% plot parameters: [vertical, horizontal] gap between subplots, margins at [bottom, top], [left, right]
height_width_gap=[0.08 0.03]; bott_top_marg =[0.05 0.05]; left_right_marg=[0.04 0.01];
params_tight_subplots={height_width_gap bott_top_marg left_right_marg};
% plot_param_settings: [fontsize_axes,fs_title,fs_legend,linewidth,params_tight_subplots,model_name]
plot_param_settings={14,22,12,4,params_tight_subplots,model_name};
% plotting stater or variables (nodes)?
state_or_node_flags={'nodes','states'};
% cutoff for minimal variation to show a variable
diff_cutoff=0.125;
figure('name','onedim parscan by param')
% fig_filename
[~,onedim_paramscan_output_cell]=fcn_onedim_parscan_plot_by_params(state_or_node_flags{2},...
		stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
		nonzero_states_inds,parscan_matrix,nodes,...
		scan_params,scan_params_up_down,... % selected parameters
		diff_cutoff,... % minimal variation for variable to be shown on plot
		plot_param_settings);

% %% SAVE
% 
% fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(diff_cutoff),'.','p')),...
% 	plot_save_folder,fig_file_type{1},'overwrite','-r200');

%% PLOT 1-dimensional scan grouped by model variables/states

% nonzero states of the model
nonzero_states_inds=find(stat_sol>0);
% sensit_cutoff: minimal value for local sensitivity or variation of model/state values
sensit_cutoff=0.1;
% parameters of plot
height_width_gap=[0.1 0.04]; bott_top_marg=[0.03 0.04]; left_right_marg=[0.07 0.02]; 
params_tight_subplots={height_width_gap bott_top_marg left_right_marg};
% plot_param_settings: [fontsize_axes,fontsize_title,params_tight_subplots,model_name]
plot_param_settings={20,20,params_tight_subplots,model_name,'colorbar'};
% plot_param_settings={12,14,[],model_name}; 
% select type of plot
plot_types={{'lineplot','heatmap'} {'nodes','states'} {'values','sensitivity'}};
plot_type_options=[1 2 1];

figure('name','1D parscan (vars)')
[resp_coeff,scan_params_sensit,scan_params_up_down_sensit,fig_filename]=...
	fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
		stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
		nonzero_states_inds,parscan_matrix,nodes,...
		scan_params,scan_params_up_down,...
		sensit_cutoff,plot_param_settings);

%%% SAVE
% fcn_save_fig(strcat(fig_filename,'_cutoff',strrep(num2str(sensit_cutoff),'.','p')),...
% 	plot_save_folder,fig_file_type{1},'overwrite','-r200');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER FITTING: define data, initial guess

% scan_params_sensit=[11 13 15 16]; scan_params_up_down_sensit={2,1,2,[1 2]};
% names of selected transition rates
[~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 

% DEFINE PARAMETER VECTOR (generate some data OR load from elsewhere)
stdev=3; meanval=0; data_param_vals=lognrnd(meanval,stdev,1,numel(predictor_names)); 
% INITIAL GUESS for parameters
init_par_vals=lognrnd(meanval,stdev,size(predictor_names)); 

% fitting NODES or STATES
var_type_flag='vars'; % 'vars' 'states'
[y_data,y_init_pred,init_error]=fcn_param_fitting_data_initguess_error(var_type_flag,...
                                   x0,stg_cell,data_param_vals,init_par_vals, stg_sorting_cell,nodes,predictor_names);
% OR TAKE 'y_data' from data

% init error:
% MSE: ((y_data-y_init_pred)'*(y_data-y_init_pred))/sum(y_data>0)
% MAE: mean(abs(y_data-fcn_statsol_values(init_par_vals)))

%%% function handles for fitting (y_data is an input!!)
[fcn_statsol_sum_sq_dev,fcn_statsol_values]=...
    fcn_handles_fitting(var_type_flag,y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);

% init_error=fcn_statsol_sum_sq_dev(init_par_vals); y_init_pred=fcn_statsol_values(init_par_vals); y_data=fcn_statsol_values(data_param_vals);

%% fitting by LSQNONLIN (local search algorithm based on trust regions)

% % [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options) 
% create fcn handles
% [~,fcn_statsol_values]=fcn_handles_fitting('states',y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);
% % feed in the function, NOT the sum of squares
% fcn_diff_statsol_data = @(p,y) y-fcn_statsol_values(p);
% y=y_data; lsqnonlin_opts=optimoptions('lsqnonlin','Display','iter');
% [bestfit_par_vals,resnorm,residual,exitflag,output]=lsqnonlin(@(p) fcn_diff_statsol_data(p,y),init_par_vals,lbnds,upbnds,lsqnonlin_opts);

% store output of lsqnonlin
lsqnonlin_opts={'iter',[],[],1e-8,[],[]}; % Display, MaxFunctionEvaluations, MaxIterations, FunctionTolerance, OptimalityTolerance, StepTolerance
lbnds=zeros(size(init_par_vals)); upbnds=1e2*ones(size(init_par_vals));
tic; [bestfit_par_vals,~,exitflag,output,history]=exastolog_lsqnonlin(var_type_flag,y_data,init_par_vals,lbnds,upbnds,...
                                                x0,stg_cell,stg_sorting_cell,nodes,predictor_names,lsqnonlin_opts); toc
% final error: fcn_statsol_sum_sq_dev(bestfit_par_vals)
%
% true values, initial guess, optimized paramset
data_init_optim_pars=[data_param_vals;init_par_vals;bestfit_par_vals];
if strcmp(var_type_flag,'states')
    data_init_optim_vars=[y_data'; y_init_pred'; fcn_statsol_values(bestfit_par_vals)']; 
    data_init_optim_vars=full(data_init_optim_vars(:,sum(data_init_optim_vars)>0));
else
    data_init_optim_vars=[y_data; y_init_pred; fcn_statsol_values(bestfit_par_vals)];
end

%% plot error & params convergence

figure('name','lsqnonlin (pars)')
%%%%% SUBPLOT1
subpl1=subplot(1,2,1); semilogy(1:size(history,1),history(:,end),'Marker','o','MarkerSize',8); xlimvals([1 size(history,1)*1.02]); set(gca,'FontSize',16)
title('SSE convergence','FontSize',24); set(subpl1,'Position',[0.05 0.08 0.44 0.86]); xlabel('iterations'); grid on
%%%%% SUBPLOT2: plot convergence of params
subpl2=subplot(1,2,2); plt_parvals=semilogy(size(history,1),data_param_vals,'Marker','o','MarkerSize',10,'LineStyle','none'); 
xlimvals([1 size(history,1)*1.02]); set(gca,'FontSize',16); 
for k=1:numel(plt_parvals); plt_parvals(k).MarkerFaceColor=plt_parvals(k).Color; end; title('parameter convergence','FontSize',24)
set(gca,'ColorOrderIndex',1); hold on; set(subpl2,'Position',[0.53 0.08 0.44 0.86]); grid on; xlabel('iterations');
semilogy(1:size(history,1),history(:,1:end-1),'LineWidth',4); legends_str=legend('NumColumns',2,'FontSize',16);
hold off;

for k=1:numel(data_param_vals)
    legends_str.String{k}=strcat('true val ',num2str(k)); legends_str.String{k+numel(data_param_vals)}=strcat('par. est. ',num2str(k)); 
end

% export_fig(strcat('colosys_plots/lsqnonlin_convergence_',num2str(numel(data_param_vals)),'params.png'),'-transparent','-nocrop','-r70')

%% plot error & variable convergence

plot_settings=[20 22];
if strcmp(var_type_flag,'vars'); sel_nodes=3:15; else; sel_nodes=[]; end
% if its states you fitted, take the transpose of ydata
% data_init_optim_vars=[y_data'; y_init_pred'; y_optim']; 
figure('name','LSQNONLIN (vars)') 
% state_var_flags={'state','var'};
error_conv=history(:,end);
fcn_plot_paramfitting(var_type_flag,data_init_optim_vars,error_conv,nodes,sel_nodes,[],[],plot_settings)
% set(gca,'xscale','log')

%%% SAVE
% fig_name=strcat('lsqnonlin_',var_type_flag,'_',num2str(numel(predictor_names)),'fittingpars');
% fcn_save_fig(fig_name,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi)

%% multidim sampling of initial conds for parameters

scan_vals=logspace(-2,2,2); n_scanvals=numel(scan_vals); n_pars=numel(predictor_names);
paramsample_table=scan_vals(fliplr(rem(floor([0:((n_scanvals^n_pars)-1)].'*n_scanvals.^(0:-1:-n_pars+1)),n_scanvals))+1);

fcn_diff_statsol_data=@(p,y) y-fcn_statsol_values(p);  y=y_data; % lsqnonlin_opts=optimoptions('lsqnonlin','Display','iter');
lsqnonlin_options=optimoptions('lsqnonlin'); paramoptim_table=zeros(size(paramsample_table) + [0 1]);

tic;
for k=1:size(paramsample_table,1)
    [bestfit_par_vals,resnorm,residual,~,~]=lsqnonlin(@(p) fcn_diff_statsol_data(p,y),paramsample_table(k,:),lbnds,upbnds,lsqnonlin_options);
    paramoptim_table(k,:)=[bestfit_par_vals resnorm];
    disp(k)
end
toc;

% 256 sets (8 params): 1133 sec ≈ 19 min -> ~4-5 sec/paramset
% 1070 sec

%% boxplot results of multistart search

figure('name','boxplot')
x_w=0.92; y_h=0.45;
% SUBPLOT1
subpl1=subplot(2,1,1); semilogy(1:size(paramoptim_table,1),paramoptim_table(:,end)/full(sum(y>0))); ylabel('MSE','FontSize',24); 
set(subpl1,'Position',[0.05 0.54 x_w y_h]); xlabel('# paramset','FontSize',24); xlimvals([0 size(paramoptim_table,1)+1]); grid on
% SUBPLOT2: boxplot of param estimates
boxplot_labels={''}; for k=1:numel(data_param_vals); boxplot_labels{k}=predictor_names{k}; end
subpl2=subplot(2,1,2); boxplot(paramoptim_table(:,1:numel(bestfit_par_vals)),'Notch','on','Labels',boxplot_labels); ylabel('fitted values','FontSize',24)
set(subpl2,'Position',[0.05 0.03 x_w 0.44]); ylim([min(min(paramoptim_table(:,1:end-1))),max(max(paramoptim_table(:,1:end-1)))])
set(gca,'YScale','log'); hold on; grid on
plot(1:numel(data_param_vals),data_param_vals,'LineStyle','none','Marker','o','MarkerFaceColor',[1 0 0],'MarkerSize',9)
legend('true value','FontSize',24)

save_flag=0; plot_save_folder='colosys_plots/';
if save_flag==1
    fig_name=strcat('lsqnonlin_',var_type_flag,'_',num2str(numel(predictor_names)),'fitpars_initconds_boxplots_parset2');
    fcn_save_fig(fig_name,plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
end

%% histograms of estim param vals vs true val

figure('name','fitpars_histogrs')
for k=1:numel(bestfit_par_vals)
    subplot(2,3,k); histogram(paramoptim_table(:,k),logspace(-4,2,10)); set(gca,'xscale','log'); % xlim([1e-2,2e2]); 
    hold on; bar(data_param_vals(k),64,10^(log10(data_param_vals(k))-1),'r','EdgeColor','none')
    title(strcat(strrep(predictor_names(k),'_','\_'),'=',num2str(round(data_param_vals(k),2)) ),'FontSize',22); ylim([0,40])
end
hold off

save_flag=0; plot_save_folder='colosys_plots/';
if save_flag==1
    fig_name=strcat('lsqnonlin_',var_type_flag,'_',num2str(numel(predictor_names)),'fitpars_histogrs');
    fcn_save_fig(fig_name,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi)
end

% param_pairs_table=nchoosek(1:size(paramoptim_table,2)-1,2);
%
% for k=1:size(param_pairs_table,1)
%   subplot(5,6,k); scatter(paramoptim_table(:,param_pairs_table(k,1)), paramoptim_table(:,param_pairs_table(k,2)))
%   set(gca,'XScale','log'); set(gca,'YScale','log'); title(strcat(num2str(param_pairs_table(k,1)),', ',num2str(param_pairs_table(k,2))))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTEGRATE cell line data

%%%%%%%%%%
% viability 

% /home/mihalykoltai/Desktop/research/models/colosys_model/data/viability/viab_df_tidy_all.csv
% viab_df_tidy_all.Properties.VariableNames
viab_df_tidy_all=readtable('../data/viability/viab_df_tidy_all.csv');
chk1i_uniqvals=unique(viab_df_tidy_all.chk1i); mk2i_uniqvals=unique(viab_df_tidy_all.mk2i); 
% we only take 0 and highest values of the inhibitors
inhib_inds=[1 1;numel(chk1i_uniqvals) 1;1 numel(mk2i_uniqvals); numel(chk1i_uniqvals) numel(mk2i_uniqvals)]; % {[0 0],[1 0],[0 1],[1 1]};
% take conditions where inhibitors are 0
cell_line_inhib_inds=cell2mat(arrayfun(@(k) find(viab_df_tidy_all.chk1i==chk1i_uniqvals(inhib_inds(k,1)) & ...
    viab_df_tidy_all.mk2i==mk2i_uniqvals(inhib_inds(k,2))),1:size(inhib_inds,1),'UniformOutput',false)');
viab_df_tidy_minmaxinhib=viab_df_tidy_all(cell_line_inhib_inds,:);

colnames=viab_df_tidy_minmaxinhib.Properties.VariableNames;
cell_line_names=unique(viab_df_tidy_minmaxinhib.cell_line); viab_per_cell_line_minmaxinhib=cell(numel(cell_line_names));
% Calculate means and stdevs grouped by cell line and conditions
viab_means_minmaxinhib=grpstats(viab_df_tidy_minmaxinhib,{'cell_line','time','mk2i','chk1i'},{'mean','std'},'DataVars','value');

% visualize cornerpoints at 1 timepoint
timeval=48; table_timeval=viab_means_minmaxinhib(viab_means_minmaxinhib.time==timeval,:);
figure('name',strcat('viab ',num2str(timeval))); ncol=7; nrow=3;
% for k=1:numel(cell_line_names)
%     indiv_cell_line=table_timeval(ismember(table_timeval.cell_line,cell_line_names{k}),:);
%     subpl_handle=subplot(nline,ncol,k);
%     heatmap(1-flipud(reshape(indiv_cell_line.mean_value,2,2)),unique(indiv_cell_line.mk2i),flipud(unique(indiv_cell_line.chk1i)),...
%         '%0.2f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,...
%         'GridLines','-','FontSize',16,'ShowAllTicks',true); set(subpl_handle,'FontSize',14); title(strrep(unique(indiv_cell_line.cell_line),'_','\_'))
%     if k>(nline-1)*ncol; xlabel('mk2'); end; if rem(k,ncol)==1; ylabel('chk1'); end
% end
ha=tight_subplot(nrow,ncol,[0.06,0.04],[0.06,0.05],[0.03,0.005]);
% [ha,pos]=tight_subplot(nrow,ncol, gap, marg_h, marg_w);
for k=1:numel(cell_line_names)
    indiv_cell_line=table_timeval(ismember(table_timeval.cell_line,cell_line_names{k}),:);
    axes(ha(k)); % subplot(nrow,ncol,k); subpl_handle=
    heatmap(1-flipud(reshape(indiv_cell_line.mean_value,2,2)),unique(indiv_cell_line.mk2i),flipud(unique(indiv_cell_line.chk1i)),...
        '%0.2f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,...
        'GridLines','-','FontSize',16,'ShowAllTicks',true); set(gca,'FontSize',14); 
    title(strrep(unique(indiv_cell_line.cell_line),'_','\_'),'FontSize',12)
    if k>(nrow-1)*ncol; xlabel('mk2'); end; if rem(k,ncol)==1; ylabel('chk1'); end
end

% export_fig(strcat('colosys_plots/celllines',num2str(timeval),'h_minmaxpoints.png'),'-transparent','-nocrop') % ,'-r70'

%% fitting model to viab data

initial_fixed_nodes={'cc','KRAS','DSB','CDC25B','g2m_trans','cell_death','CHEK1','MAPKAPK2'}; 
% [chk1i,mk2i]={[0,0],[1,0],[0,1],[1,1]}
inhib_combs=1-[0,0;1,0;0,1;1,1]; initial_fixed_nodes_vals=[1,1,1,0,0,0];
initial_fixed_nodes_vals_matr=[repmat(initial_fixed_nodes_vals,4,1),inhib_combs]; stationary_node_vals_perturbs=zeros(4,numel(nodes)); 
% init_node_vals_perturbs=stationary_node_vals_perturbs;
% stg_cell=fcn_build_stg_cell(truth_table_filename,nodes);

% stationary_node_vals_perturbs(:,ismember(nodes,{'DSB','CHEK1','MAPKAPK2','g2m_trans','cell_death'}))
% [~,d]=ismember(initial_fixed_nodes,nodes); array2table(stationary_node_vals_perturbs(:,d),'VariableNames',initial_fixed_nodes)

viability_data=1-indiv_cell_line.mean_value';
% node (or state) to fit
fit_vars=10; perturb_nodes={'CHEK1','MAPKAPK2'};
distrib_types={'random','uniform'};distr_typ_sel=2;
[x,resnorm,exitflag,output,history]=exastolog_lsqnonlin_simult_fitting(var_type_flag,viability_data,init_par_vals,lbnds,upbnds,...
                                                stg_cell,nodes,perturb_nodes,predictor_names,... % predictor_names,
                                                fit_vars,...
                                                distrib_types,distr_typ_sel,dom_prob,...
                                                initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,lsqnonlin_opts);

% init_par_vals=[5,0.8287,0.0169]

%% multidim parscan

predictor_names={'d_DSB','u_g2m_trans','u_CDC25B'}; % {'d_DSB','u_CHEK1','u_g2m_trans','d_HRR2','u_MAPKAPK2','u_CDC25B'}; 
distrib_types={'random','uniform'};distr_typ_sel=2;
scan_vals=logspace(-2,2,5); n_scanvals=numel(scan_vals); n_pars=numel(predictor_names);
paramsample_table=fcn_paramsample_table(scan_vals,n_pars);
% obj fcn
fcn_perturb_exps_optim=@(x)objfcn_perturb_exps(var_type_flag,fit_vars,stg_cell,nodes,perturb_nodes,...
                distrib_types,distr_typ_sel,dom_prob,initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,...
                fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x));
output_perturb=zeros(size(paramsample_table,1),5);
% scan
for k=1:size(paramsample_table,1)
    simul_y=fcn_perturb_exps_optim(paramsample_table(k,:)); 
    output_perturb(k,:)=[simul_y,sum((viability_data-simul_y).^2)];
    if rem(k,10)==0; disp(k); end
end

%% plot error as fcn of par values

figure('name','sse-multiscan')

subpl1=subplot(2,3,1); plot(output_perturb(:,end)/3,'Marker','o'); xlabel('parameter sets','FontSize',12); ylabel('mean squared error','FontSize',12); 
set(subpl1,'Position',[0.13 0.58 0.78 0.34])
for k=1:size(paramsample_table,2)
subplot(2,3,k+3)
    scatter(paramsample_table(:,k),output_perturb(:,end)); set(gca,'xscale','log')
    title(predictor_names(k),'Interpreter','latex','FontSize',14)
    hold on; sampl_means=cell2mat(arrayfun(@(x) mean(output_perturb(paramsample_table(:,k)==x,end)),unique(paramsample_table(:,k)),'UniformOutput',false));
    sampl_std=cell2mat(arrayfun(@(x) std(output_perturb(paramsample_table(:,k)==x,end)),unique(paramsample_table(:,k)),'UniformOutput',false));
    plot(unique(paramsample_table(:,k)),sampl_means,'Marker','o','Color','red','MarkerFaceColor','red'); 
    errorbar(unique(paramsample_table(:,k)),sampl_means,sampl_std,'CapSize',12,'Color','red','LineStyle','none','LineWidth',2)
    hold off; if rem(k,3)==1; ylabel('SSE','FontSize',12); end
end

% fig_name=strcat('error_parscan_',strcat(num2str(numel(predictor_names)),'pars')); 
% fcn_save_fig(fig_name,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi)

%% regression: error ~ paramsample_table

% on error
% linregr_err=fitlm(paramsample_table,output_perturb(:,end));
% plot(1:size(paramsample_table,1),[output_perturb(:,end),predict(linregr_err)]); xlim([1,size(paramsample_table,1)+2]); legend({'boolean model error','prediction'})

% regress on model output
% n_var=3;
figure('name','allconderrs_regr')
for n_var=2:4
linregr_output=fitlm(paramsample_table,output_perturb(:,n_var)); linregr_pred=predict(linregr_output);
% nonlinear fit
fcn_modred=@(x,p) sum(( p(:,4) - (p(:,3)./(p(:,3)+x(3))).*( p(:,2)./(p(:,2)+x(2)) ).*( 1./(p(:,1).^x(4)+x(1) ) ) ).^2);
p=[paramsample_table,output_perturb(:,n_var)];
[X,FVAL,EXITFLAG] =fminsearch(@(x) fcn_modred(x,p),[10,10,0.1,2]);
x1=paramsample_table(:,1); x2=paramsample_table(:,2); x3=paramsample_table(:,3); 
c1=X(1); c2=X(2); c3=X(3);  c4=X(4); nonlin_predvals=(x3./(x3+c3)).*(x2./(x2+c2)).*(1./(x1.^c4+c1));
nonlin_pred_legend=strcat('algebr fit: (p_3/(p_3+',num2str(round(c3,2)),'))*(p_2/(p_2+',...
    num2str(round(c2,2)),'))*(1/(p_1^{',num2str(round(c4,2)),'}+',num2str(round(c1,2)),'))');
% PLOT
% 
subpl1=subplot(4,1,n_var-1); plot(output_perturb(:,n_var),'Marker','o','MarkerFaceColor','blue'); hold on; plot(nonlin_predvals)
legend({strcat('model output (cell\_death) condition ',num2str(n_var)),nonlin_pred_legend},'FontSize',16)
xlimvals([1,size(paramsample_table,1)+2]); xticks_val=sort([5:5:size(paramsample_table,1),(0:5:120)+1]); xtick_label_vals=cell(1,size(paramsample_table,1)); 
for k = 1:numel(xticks_val); xtick_label_vals{xticks_val(k)}=xticks_val(k); end
xticks(1:125); xtickangle(90); xticklabels(xtick_label_vals); % set(subpl1,'Position',[0.05 0.54 0.9 0.44]); grid on
end
%%%
subpl2=subplot(4,1,4); heatmap(log10(paramsample_table)','',predictor_names,'','GridLines','-'); colorbar; xlabel('log10(parval)'); 
% set(subpl2,'Position',[0.05 0.03 0.9 0.44])

% fcn_save_fig('model_output_lin_regr',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
% fcn_save_fig('model_output_algebr_fit',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%% without perturbs

scan_vals=logspace(-2,2,4); n_scanvals=numel(scan_vals); n_pars=numel(predictor_names);
paramsample_table=fcn_paramsample_table(scan_vals,n_pars);
% [~,~,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes);
[stat_sol_paramsample_table,stat_sol_states_paramsample_table]=fcn_calc_paramsample_table(paramsample_table,...
                                                                scan_params_sensit, scan_params_up_down_sensit,...
                                                                transition_rates_table,stg_cell,x0,100);


%% regression: error/model output ~ paramsample_table

% on error
% linregr_err=fitlm(paramsample_table,output_perturb(:,end));
% plot(1:size(paramsample_table,1),[output_perturb(:,end),predict(linregr_err)]); xlim([1,size(paramsample_table,1)+2]); legend({'boolean model error','prediction'})

% regress on model output
% which var changes the most?
std_per_node=arrayfun(@(X) std(stat_sol_paramsample_table(:,X)),1:size(stat_sol_paramsample_table,2));
[~,n_var]=max(std_per_node); % find(ismember(nodes,'cell_death'));
output_perturb=stat_sol_paramsample_table;

linregr_output=fitlm(paramsample_table,output_perturb(:,n_var)); linregr_pred=predict(linregr_output);

% nonlinear fit
fcn_modred=@(fitpars,inputpars) sum( ( inputpars(:,end) - ...
   ( ...
     ( inputpars(:,2).^fitpars(1) + inputpars(:,3).^fitpars(2) + inputpars(:,4).^fitpars(3) + inputpars(:,5).^fitpars(4) + inputpars(:,6).^fitpars(5) )./...
    ( inputpars(:,2).^fitpars(6)+ inputpars(:,3).^fitpars(7) + inputpars(:,4).^fitpars(8) + inputpars(:,5).^fitpars(9) + inputpars(:,6).^fitpars(10) )  ...
    ) ).^2);
% simpler model 
% hyperboles
% 

inputpars=[paramsample_table,output_perturb(:,n_var)];

% with FMINSEARCH
fminsrch_optimset=optimset; fminsrch_optimset.Display='final'; fminsrch_optimset.MaxFunEvals=1e4; fminsrch_optimset.TolFun=1e-6; fminsrch_optimset.MaxIter=1e4;
fminsrch_optimset.TolX=1e-6; 
initconds_lhs=lhsdesign(100,10);
for k=1:size(initconds_lhs,1)
	x0_fminsrch=abs(randn(1,10)); x0_fminsrch=x0_fminsrch/sum(x0_fminsrch);
    % x0_fminsrch=initconds_lhs(k,:);
    [X,FVAL,EXITFLAG]=fminsearch(@(fitpars) fcn_modred(fitpars,inputpars),x0_fminsrch,fminsrch_optimset);
    rmse_nonlin=sqrt(FVAL/size(paramsample_table,1));
    initparvals_optparvals_error_fminsrch(k,:)=[x0_fminsrch,X,FVAL];
    disp(k)
end 

% best param set
[~,b]=min(initparvals_optparvals_error_fminsrch(:,end)); fitpars=initparvals_optparvals_error_fminsrch(b,(1:numel(X))+numel(X)); 
rmse_nonlin=sqrt(initparvals_optparvals_error_fminsrch(b,end)/size(paramsample_table,1));

% with LSQNONLIN
fcn_modred_lsqnonlin=@(fitpars) inputpars(:,end) - ...
 (  fitpars(1)*( inputpars(:,2).^fitpars(10))./(inputpars(:,2).^fitpars(10) + fitpars(2)) + ...
    fitpars(3)*1./(inputpars(:,3).^fitpars(11) + fitpars(4)) + ...
    fitpars(5)*1./(inputpars(:,4) + fitpars(6)) + ...
    fitpars(7)*(inputpars(:,5).^fitpars(12))./(inputpars(:,5).^fitpars(12) + fitpars(8)) + ...
    fitpars(8)*1./(inputpars(:,6) + fitpars(9)) );
n_pars=12; lsqnonlin_opts=optimoptions('lsqnonlin'); lsqnonlin_opts.Display='final';
for k=1:size(initconds_lhs,1)
    disp(k)
    [X,resnorm,residual,exitflag,output]=lsqnonlin(fcn_modred_lsqnonlin,randn(1,n_pars),[],[],lsqnonlin_opts); % zeros(1,n_pars),10*ones(1,n_pars)
	initparvals_optparvals_error_lsqnonlin(k,:)=[X,resnorm];
    resnorm
end

% best param set
[~,b]=min(initparvals_optparvals_error_lsqnonlin(:,end)); fitpars=initparvals_optparvals_error_lsqnonlin(b,(1:numel(X))); 
rmse_nonlin=sqrt(initparvals_optparvals_error_lsqnonlin(b,end)/size(paramsample_table,1));

nonlin_predvals=fitpars(1)*( inputpars(:,2).^fitpars(10))./(inputpars(:,2).^fitpars(10) + fitpars(2)) + ...
    fitpars(3)*1./(inputpars(:,3).^fitpars(11) + fitpars(4)) + ...
    fitpars(5)*1./(inputpars(:,4) + fitpars(6)) + ...
    fitpars(7)*(inputpars(:,5).^fitpars(12))./(inputpars(:,5).^fitpars(12) + fitpars(8)) + ...
    fitpars(8)*1./(inputpars(:,6) + fitpars(9));
% check if correct
% sum((output_perturb(:,n_var) - nonlin_predvals).^2)
R_sq_nonlin=1-sum((output_perturb(:,n_var)-nonlin_predvals).^2)/sum( (output_perturb(:,n_var)-mean(output_perturb(:,n_var)) ).^2);
R_sq_lin=linregr_output.Rsquared.Ordinary;

figure('name','parfit')
xlimvals=1:1000; subpl1=subplot(2,1,1);
plot(output_perturb(:,n_var),'Marker','o','LineStyle',':','MarkerSize',4); hold on; plot(nonlin_predvals) % plot(linregr_pred(xlimvals)); 
legend({'model output','nonlin regr'},'FontSize',18); % set(gca,'Position',[0.044,0.316,0.94,0.67]); hold on
title(strcat('nonlinear alg eq R^2=',num2str(round(R_sq_nonlin,2)),', RMSE=',num2str(round(rmse_nonlin,3)),' (lin regr R^2=',num2str(round(R_sq_lin,2)),')'),'FontSize',18)
xlim(xlimvals([1,end])); % heatmap of param tables
subpl2=subplot(2,1,2); 
heatmap(log10(paramsample_table(xlimvals,:))','',predictor_names,'',...
    'Colormap','redblue','MinColorValue',min(log10(paramsample_table(:))),'MaxColorValue',max(log10(paramsample_table(:))),'GridLines','-'); 
xlabel('log10(parval)'); set(subpl2,'Position',[0.04, 0.0316,0.94,0.241])
set(subpl1,'Position',[0.044,0.316,0.94,0.64]);

% fcn_save_fig('model_output_lin_regr',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)
% fcn_save_fig('model_output_algebr_fit',plot_save_folder,fig_file_type{3},'overwrite',resolution_dpi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WES (phosphoprotein)

% /home/mihalykoltai/Desktop/research/models/colosys_model/data/wes_facs/wes_facs_data_tidy_scaled_vals.csv
% /home/mihalykoltai/Desktop/research/models/colosys_model/data/wes_facs/summary_WES_data_and_cell_line_data.csv
% /home/mihalykoltai/Desktop/research/models/colosys_model/data/wes_facs/wes_facs_data.csv
% /home/mihalykoltai/Desktop/research/models/colosys_model/data/wes_facs/wes_facs_data_tidy.csv

% as a table
wesdata_filename='../data/wes_facs/wes_facs_data_tidy_scaled_vals.csv'; readtable_opts=detectImportOptions(wesdata_filename);
value_cols=cell2mat(cellfun(@(x) contains(x,'value'),readtable_opts.VariableNames,'UniformOutput',false));
wes_facs_data_tidy_scaled_vals=readtable(wesdata_filename,readtable_opts);
wes_exp_inds=ismember(wes_facs_data_tidy_scaled_vals.expertype,'wes');
sel_vars={'cell_line','trunc_variable','chek1i','mk2i','value','log2_scaled_value',...
           'scaled_value','scaled_value_no_etop','scaled_value_no_etop_across_samples'};
% sel_vars_inds=ismember(wes_facs_data_tidy_scaled_vals.Properties.VariableNames, sel_vars);
[~,sel_vars_inds]=ismember(sel_vars,wes_facs_data_tidy_scaled_vals.Properties.VariableNames);
wes_table=wes_facs_data_tidy_scaled_vals(wes_exp_inds,sel_vars);
% as a matrix
wes_facs_data_tidy_scaled_vals_matrix=readmatrix(wesdata_filename); wes_matrix=wes_facs_data_tidy_scaled_vals_matrix(wes_exp_inds,sel_vars_inds);
% {'chek1i','mk2i','etoposide','treatment','cell_line','expertype','time','variable',
% 'value','log2_scaled_value','scaled_value','scaled_value_no_etop','scaled_value_no_etop_across_samples','trunc_variable'}

wes_table_vars=unique(wes_table.trunc_variable); wes_table_cell_lines=unique(wes_table.cell_line);

% visualize on heatmap 4 values for the variables per cell line
n_ln=numel(wes_table_cell_lines); n_col=numel(wes_table_vars);
figure('name','wes')
for k_cell_line=1:n_ln
for k_var=1:n_col
    indiv_vals=wes_table(ismember(wes_table.trunc_variable,wes_table_vars(k_var)) & ismember(wes_table.cell_line,wes_table_cell_lines(k_cell_line)),:);
    chk1ivals=unique(indiv_vals.chek1i); mk2ivals=unique(indiv_vals.mk2i);    
    pltcntr=(k_cell_line-1)*n_col+k_var; subpl_handle=subplot(n_ln,n_col,pltcntr);
    indiv_vals_matr=reshape(indiv_vals.scaled_value_no_etop,numel(chk1ivals),numel(mk2ivals));
    % [1 end]
    heatmap(flipud(indiv_vals_matr([1 end],:)),mk2ivals,flipud(chk1ivals([1 end])),...
        '%0.2f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,...
        'GridLines','-','FontSize',16,'ShowAllTicks',true); set(subpl_handle,'FontSize',14); title(strrep(wes_table_vars(k_var),'_','\_'))
    if rem(k_var,n_col)==1; ylabel({wes_table_cell_lines{k_cell_line},'[chek1]'}); end
    if pltcntr>n_col*(n_ln-1); xlabel('[mk2]'); end
end
end

save_flag=0; 
if save_flag==1
filetype='.png';
if size(indiv_vals_matr,1)==2
    export_fig(strcat('colosys_plots/wes_minmaxpoints',filetype),'-transparent','-nocrop') % ,'-r70'
else
    export_fig(strcat('colosys_plots/wes_allpoints',filetype),'-transparent','-nocrop') % ,'-r70'
end
end

%% Cytof (phosphoprotein)

% 4 cell lines, there are two files RAW and LOG2 values
% ~/Desktop/research/models/colosys_model/data/cytof/All means_DMSO_alllog2.csv
% ~/Desktop/research/models/colosys_model/data/cytof/All means_DMSO_allraw.csv

% protein names:
protein_names={'pChk2','cCasp3','cPARP','pan_Akt','Axin2','PTK7','pH2A.X','pChk1','pSmad1/8','YAP','pMEK1/2','pAkt','pSmad2/3',...
'PROM1','pNF-κB','p-p38','EphB2','pCDC25c','LGR5','Ki-67','IκBα','Vimentin','CD44','EpCAM','pAkt_T308','CD24','p4e-BP1','pERK1/2','p-p53',...
'act Notch-1','total_ERK','pS6'};

means_DMSO_allraw_filename='../data/cytof/All means_DMSO_allraw.csv'; readtable_opts=detectImportOptions(means_DMSO_allraw_filename);
value_cols=cell2mat(cellfun(@(x) contains(x,'value'),readtable_opts.VariableNames,'UniformOutput',false));
means_allraw=readtable(means_DMSO_allraw_filename,readtable_opts);
cytof_means_all_log2=readtable('../data/cytof/All means_DMSO_alllog2.csv',readtable_opts);
cell_line_names=unique(cytof_means_all_log2.cell_line);
% heatmap
figure('name','cytof'); n_col_vals=5:size(cytof_means_all_log2,2);
% select variables that show some minimal variation
% resp_vars=max(table2array(means_all_log2(:,n_col_valstart:end))) - min(table2array(means_all_log2(:,n_col_valstart:end)))>0.15;
grprows=arrayfun(@(x) (1:numel(cell_line_names))+(x-1)*numel(cell_line_names),1:numel(cell_line_names),'UniformOutput',false); 
grptbls=cellfun(@(x) cytof_means_all_log2(x(:),n_col_vals),grprows,'UniformOutput',false);
grptbls_minmax=cell2mat(arrayfun(@(x) max(table2array(grptbls{x}))-min(table2array(grptbls{x})),1:numel(grptbls),'UniformOutput',false)');
thrsh_val=0.15; selvars=sum(grptbls_minmax>thrsh_val)>0;
for k=1:numel(cell_line_names)
    if k==1; xlab_input={'Chk1i','Mk2i'}; else; xlab_input=''; end
    subpl1=subplot(k,1,1); set(subpl1,'Position',[0.01,0.14+(k-1)*0.22,0.08,0.18]); heatmap([0,0;1,0;0,1;1,1],xlab_input,'',...
        '%0.0f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,'GridLines','-','FontSize',16,'ShowAllTicks',true); 
    set(subpl1,'FontSize',18); subpl2=subplot(k,2,1); set(subpl2,'Position',[0.12,0.14+(k-1)*0.22,0.87,0.18])
    inds=ismember(cytof_means_all_log2.cell_line,cell_line_names(k));
    if k==1; xlabs=protein_names(selvars); else; xlabs=''; end % means_all_log2.Properties.VariableNames(n_col_vals(selvars))
    conds=cellfun(@(x) strsplit(x,'_'),cytof_means_all_log2.cell_line_cond(inds),'UniformOutput',false); % arrayfun(@(x) conds{x}(end), 1:numel(conds))
    heatmap(table2array(cytof_means_all_log2(inds,n_col_vals(selvars))),xlabs,'',...
        '%0.2f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,'GridLines','-','FontSize',10,'ShowAllTicks',true); 
    set(subpl2,'FontSize',14); ylabel(cell_line_names(k),'FontSize',18)
end

filetype='.pdf'; save_flag=0; 
if save_flag==1; export_fig(strcat('colosys_plots/cytof_log2fc_thrshld',num2str(thrsh_val),filetype),'-transparent','-nocrop'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fminsearch (downhill simplex method)
% 
% % [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(...);
% % f = @(x,scan_vals=logspace(0,1,5); % sort(data_param_vals)+normrnd(0,0.1,1,5); % logspace(-0.5,1,n_scanvals);
% % fminsearch options
% fminsearch_opts={'iter',2000,[]}; % 2000,2000 'Display','iter', 'MaxFunEvals', MaxIter
% % RUN DOWNHILL SIMPLEX search
% tic; [bestfit_par_vals,fval,exitflag,output,history]=exastolog_fminsearch(var_type_flag,y_data,init_par_vals,...
%                                                         x0,stg_cell,stg_sorting_cell,nodes,predictor_names,fminsearch_opts); toc
% % DATA, INIT, OPTIM parameter values
% data_init_optim_pars=[data_param_vals;init_par_vals;bestfit_par_vals];
% data_init_optim_vars=[y_data'; y_init_pred'; fcn_statsol_values(bestfit_par_vals)']; 
% data_init_optim_vars=full(data_init_optim_vars(:,sum(data_init_optim_vars)>0));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Multidimensional parameter scanning with Latin Hypercube Sampling
% 
% sampling_types={'lognorm','linear','logunif'}; sampling_type=sampling_types{3};
% % par_min_mean: minimum or (if lognormal) mean of distribution. 
% %  				Can be a scalar or vector (if different values for different parameters)
% % max_stdev: maximum or in case of lognormal the mean of distribution. Scalar or vector
% %
% % for 'lognorm','logunif' provide LOG10 value of desired mean/min & stdev/max (-2 is a mean of 0.01)
% par_min_mean=-2; % repmat(1.5,1,numel(cell2mat(scan_params_up_down_sensit(:)))); par_min_mean(4)=3; 
% max_stdev=2; 	 % repmat(0.5,1,numel(cell2mat(scan_params_up_down_sensit(:))));
% % number of param sets
% lhs_scan_dim=1000;
% 
% % RUN the LHS
% [all_par_vals_lhs,stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan]=... % outputs
%     fcn_multidim_parscan_latinhypcube(par_min_mean,max_stdev,sampling_type,lhs_scan_dim, ...
%                           scan_params_sensit,scan_params_up_down_sensit, ...
%                           transition_rates_table,stg_cell,x0,nodes);
%                       
% %% Visualize LHS by scatter plots
% 
% % which variable to plot?
% var_ind=10;
% % STATES or NODES? <scan_values>: values to be plotted
% % model variables: stat_sol_nodes_lhs_parscan; states: stat_sol_states_lhs_parscan
% scan_values=stat_sol_nodes_lhs_parscan; % stat_sol_states_lhs_parscan
% 
% sampling_type=sampling_types{3}; % sampling_types={'lognorm','linear','logunif'};
% % file_name_prefix=strcat('LHS_parscan_scatterplot_trend_',nodes{var_ind}); 
% file_name_prefix=strcat('LHS_parscan_scatterplot_trend_state',num2str(var_ind));
% % param_settings: [number_bins_for_mean,trendline_width,axes_fontsize,index nonzero states]
% param_settings = [50 6 24 size(stat_sol_states_lhs_parscan)];
% 
% figure('name',num2str(var_ind))
% fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,...
%         scan_params_sensit,scan_params_up_down_sensit,nodes,sampling_type,param_settings)
% 
% % %% SAVE
% % resolution_dpi='-r200'; 
% % fcn_save_fig(file_name_prefix,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi);
% 
% %% Regression of variables by transition rates
% 
% plot_type_flag={'par_var','heatmap','r_sq'}; % {'par_var','heatmap'/'lineplot','r_sq'/'slope'}
% sel_nodes=[3:7,9:15]; % if empty --> all
% % plot_settings=[fontsize,fontsize,maximum value for heatmap colors], 
% % if plot_settings(3)=NaN, then max color automatically selected
% plot_settings=[20 20 0.29]; 
% % if regression type is 'linlog', then the fit is y = a + b*log10(x)
% regr_types={'log','linear'}; % log recommended if parameter values log-uniformly distributed in sampling
% figure('name',strjoin(plot_type_flag))
% scan_values=stat_sol_nodes_lhs_parscan; % or:  stat_sol_states_lhs_parscan
% [r_squared,slope_intercept]=fcn_multidim_parscan_parvarcorrs(plot_type_flag,...
% 		all_par_vals_lhs,scan_values,...
% 		nodes,sel_nodes,... % which nodes
% 		scan_params_sensit,scan_params_up_down_sensit, ... % same params as in LHS!
% 		regr_types{1},plot_settings);
% 
% % %% SAVE plot
% % fig_prefix=strjoin(plot_type_flag,'_'); resolution_dpi='-r350'; 
% % fcn_save_fig(fig_prefix,plot_save_folder,fig_file_type{1},'overwrite',resolution_dpi)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% error as a fcn of par vals
% 
% [sensit_params_table,sensit_params_seq,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% 
% parscan_min_max=[1e-2,1e3]; n_steps=20; sampling_types={'log','linear'};
% 
% % matrix of parameter values: 1D sweeps
% parscan_matrix=fcn_onedim_parscan_generate_matrix(scan_params_sensit,scan_params_up_down_sensit,...
% 		nodes,sampling_types{1},parscan_min_max,n_steps);
% 
% % calculation
% tic;
% [stationary_state_vals_onedimscan,stationary_node_vals_onedimscan,stationary_state_inds_scan]=...
%     fcn_onedim_parscan_calc(stg_cell,transition_rates_table,x0,nodes,parscan_matrix,scan_params_sensit,scan_params_up_down_sensit);
% toc
% 
% % sweeping in 1 parameter:
% % semilogx(parscan_matrix(:,1),squeeze(stationary_state_vals_onedimscan(1,:,:)),'LineWidth',3)
% % data: y_data'
% % error when sweeping in 1 parameter:
% sse=zeros(fliplr(size(parscan_matrix)));
% figure
% for k=1:size(parscan_matrix,2)
%     subplot(2,3,k)
%     sse(k,:)=sum((y_data'-squeeze(stationary_state_vals_onedimscan(k,:,:))).^2,2);
%     semilogx(parscan_matrix(:,k),sse(k,:),'LineWidth',3)
%     title(strrep(predictor_names(k),'_','\_')); set(gca,'FontSize',14); grid on;
% end
% subplot(2,3,6); semilogx(parscan_matrix(:,k),sum(sse)); title('all params'); set(gca,'FontSize',14); grid on;

%% multidimensional param sweep with regular grids
% 
% % if we want to scan in p parameters, with n values, this is a n^p point
% % param scan, eg. 5 params 3 values --> 3^5=243 sampling points
% 
% scan_vals=logspace(0,1,5); % sort(data_param_vals)+normrnd(0,0.1,1,5); % logspace(-0.5,1,n_scanvals);
% n_scanvals=5; n_pars=numel(scan_params_sensit);
% paramsample_table=scan_vals(fliplr(rem(floor([0:((n_scanvals^n_pars)-1)].'*n_scanvals.^(0:-1:-n_pars+1)),n_scanvals))+1);
% % [size(paramsample_table),n_scanvals^n_pars]
% 
% % scan_params_sensit,scan_params_up_down_sensit
% 
% % transition rates to scan in
% disp_var=5; % show at every n% the progress
% [stat_sol_paramsample_table,stat_sol_states_paramsample_table]=...
%     fcn_calc_paramsample_table(paramsample_table,scan_params_sensit,scan_params_up_down_sensit,...\
%     transition_rates_table,stg_cell,x0,disp_var);

%% plot error
% sse_multidim_scan=sum((y_data'-stat_sol_states_paramsample_table).^2,2);
% 
% plt_thr=1e-3;
% figure('name','parsets by SSE')
% semilogy(1:sum(sse_multidim_scan<plt_thr),sse_multidim_scan(sse_multidim_scan<plt_thr),'Marker','o','MarkerFaceColor','red'); 
% xlim([0 sum(sse_multidim_scan<plt_thr)+1])
% text(1:sum(sse_multidim_scan<plt_thr),sse_multidim_scan(sse_multidim_scan<plt_thr),num2str(find(sse_multidim_scan<plt_thr))); grid on
% 
% % xlim([1,size(stat_sol_states_paramsample_table,1)])
% [a,best_fit_parset]=min(sse_multidim_scan); % best fit
% [~,closest_parset]=min(sum((data_param_vals-paramsample_table).^2,2)); % distance from actual params
% [data_param_vals;paramsample_table([best_fit_parset,closest_parset],:)] % actual, best fit, closest to params
% 
% %% recalc with 1 parameter set
% 
% % [sensit_params_table,sensit_params_seq,predictor_names]=fcn_get_trans_rates_tbl_inds(scan_params_sensit,scan_params_up_down_sensit,nodes); 
% 
% % parsets
% % [data_param_vals;paramsample_table([best_fit_parset,closest_parset],:)] % actual, best fit, closest to params
% 
% % transition_rates_table=fcn_trans_rates_table(nodes,distr_type{1},meanval,sd_val,chosen_rates,chosen_rates_vals);
% sel_par=[2]; parscan_size=120;
% figure('name',strcat('sampling best parset', num2str(sel_par)))
% parind=closest_parset;
% for j=1:numel(sel_par)
% for k=1:parscan_size
%     
% z_vals(k)=paramsample_table(parind,sel_par(j))+(k-parscan_size/2)*0.1; if z_vals(k)<0; z_vals(k)=1e-5; end
% p_mod=paramsample_table(parind,:); p_mod(sel_par(j))=z_vals(k);
% %
% transition_rates_table(scan_params_sensit)=p_mod; % [1.1167,7.2863,3.4549,4.3764,3.7173];
% %%%%
% tic; [A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,''); toc
% [stat_sol_blocks,~,~]=split_calc_inverse(A_sparse, stg_sorting_cell, transition_rates_table, x0);
% parscan_diffs(k,:)=full(stat_sol_blocks(stat_sol_blocks>0))'; %diff([y_data';full(stat_sol_blocks(stat_sol_blocks>0))']);
% end
% 
% subplot(2,2,1+(j-1)*2); plot(z_vals,parscan_diffs,'LineWidth',3); title(strcat('deviation by states, par',num2str(sel_par(j)))); legend
% subplot(2,2,2+(j-1)*2); semilogy(z_vals,sum(parscan_diffs.^2,2),'Marker','o','MarkerFaceColor','red','LineWidth',3); title('SSE')
% end
% 
% %% error of single example
% 
% transition_rates_table(sensit_params_seq)=paramsample_table(closest_parset,:); % [1.1167,7.2863,3.4549,4.3764,3.7173];
% %%%%
% tic; [A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,''); toc
% [stat_sol_blocks,~,~]=split_calc_inverse(A_sparse, stg_sorting_cell, transition_rates_table, x0);
% d=diff([y_data';full(stat_sol_blocks(stat_sol_blocks>0))']);
% d*d'