function [figure_name,output_cell]=fcn_onedim_parscan_plot_by_params(state_or_node_flag,...
                                        stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,...
                                        nonzero_states_inds,parscan_matrix,nodes,scan_params,scan_params_up_down,...
                                        diff_cutoff,plot_param_settings)

% plot parameters
title_fontsize=plot_param_settings{2}; legend_fontsize=plot_param_settings{3}; linewidth_val=plot_param_settings{4};
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
[~,~,param_names]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
if strcmp(state_or_node_flag,'states')
    scan_variable=stationary_state_vals_onedimscan;
elseif strcmp(state_or_node_flag,'nodes')
    scan_variable=stationary_node_vals_onedimscan;
else
    error('<state_or_node_flag> has to be ''nodes'' or ''states''.')
end
sensit_params_table=arrayfun(@(x) (max(scan_variable(:,:,x),[],2)-min(scan_variable(:,:,x),[],2))'>diff_cutoff,...
    1:size(scan_variable,3),'un',0);
% sensit_params_table=vertcat(sensit_params_table{:});
sensit_params = find(sum(vertcat(sensit_params_table{:}))>0); 
% sensit_states=nonzero_states_inds(sum(vertcat(sensit_params_table{:}),2)>0)';
% variables sensitive to some parameter
sensit_vars_all=find(sum(vertcat(sensit_params_table{:}),2)>0)';
all_colors=distinguishable_colors(numel(sensit_vars_all));
nrow=round(sqrt(numel(sensit_params))); ncol=nrow; if nrow*ncol<numel(sensit_params); ncol=ncol+1; end

tight_subplot_params=plot_param_settings{5};
if ~isempty(tight_subplot_params)
    
    [ha,~]=tight_subplot(nrow,ncol,tight_subplot_params{1},tight_subplot_params{2},tight_subplot_params{3});
end
output_cell=cell(numel(sensit_params),1);

for k=1:size(scan_variable,1)
    data_plot=squeeze(scan_variable(k,:,:));
    sensit_vars_indiv_param=find(max(data_plot)-min(data_plot)>diff_cutoff);
    if k==1; counter=0; end
    if ~isempty(sensit_vars_indiv_param)
      counter=counter+1;
      if ~isempty(tight_subplot_params)
        axes(ha(counter)); 
      else
          subplot(nrow,ncol,counter)
      end
      
      color_vals = all_colors(ismember(sensit_vars_all,sensit_vars_indiv_param),:);
      for j=1:numel(sensit_vars_indiv_param)
        semilogx(parscan_matrix(:,k),data_plot(:,sensit_vars_indiv_param(j)),'LineWidth',linewidth_val,'Color',color_vals(j,:)); hold on;
      end
      % legends for states
   if strcmp(state_or_node_flag,'states')
    states_legend_index=cellfun(@num2str,num2cell(nonzero_states_inds(sensit_vars_indiv_param)),'un',0)'; 
    states_legend_binary=num2str(truth_table_inputs(nonzero_states_inds(sensit_vars_indiv_param),:));
    states_legend_binary=arrayfun(@(x) strrep(states_legend_binary(x,:),' ',''), 1:size(states_legend_binary,1), 'un', 0 );
    legend_combined_strings=arrayfun(@(x) strcat(states_legend_index{x},', [',states_legend_binary{x},']'),1:numel(states_legend_index),'un',0);
    legend(states_legend_index,'FontSize',legend_fontsize,'Interpreter','none'); % legend_combined_strings
    output_cell{counter} = {param_names(k) nonzero_states_inds(sensit_vars_indiv_param)};
   elseif strcmp(state_or_node_flag,'nodes')
          legend(nodes(sensit_vars_indiv_param),'FontSize',legend_fontsize,'Interpreter','none');
          output_cell{counter} = {param_names(k) nodes(sensit_vars_indiv_param)};
   end
      title(param_names(k),'Interpreter','none','FontSize',title_fontsize)
      % output
    end
end

model_name=plot_param_settings{6};
figure_name = strcat('onedim_parscan_lineplot_',state_or_node_flag,'_',model_name,'_by_params');