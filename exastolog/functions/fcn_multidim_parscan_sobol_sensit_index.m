function tau_i=fcn_multidim_parscan_sobol_sensit_index(sobol_sensit_index,var_type,all_par_vals_lhs,...
                                            stat_sol_nodes_lhs_parscan,stat_sol_states_lhs_parscan,...
                                            sample_size,...
                                            sequential_indices_lhs,...
                                            scan_params_sobol,scan_params_up_down_sobol,...
                                            stg_cell,transition_rates_table,x0,nodes,...
                                            sel_nodes,plot_settings,disp_freq)

par_ind_table=[repelem(scan_params_sobol, cellfun(@(x) numel(x),scan_params_up_down_sobol))', horzcat(scan_params_up_down_sobol{:})'];
prefix={'u_','d_'}; 
predictor_names=arrayfun(@(x) strcat(prefix(par_ind_table(x,2)), nodes(par_ind_table(x,1))), 1:size(par_ind_table,1),'un',0);
predictor_names=vertcat(predictor_names{:})'; 

[~,sequential_indices_sobol,~] = fcn_get_trans_rates_tbl_inds(scan_params_sobol,scan_params_up_down_sobol,nodes);

% if Sobol index calculated not for all params of original LHS scan, but only a subset
if numel(sequential_indices_sobol)~=numel(sequential_indices_lhs)
    all_par_vals_lhs_subset=all_par_vals_lhs(:,ismember(sequential_indices_lhs,sequential_indices_sobol));
else
    all_par_vals_lhs_subset=all_par_vals_lhs;
end

if strcmp(var_type,'node')
    scan_values=stat_sol_nodes_lhs_parscan;
    if isempty(sel_nodes)
        sel_nodes=1:size(scan_values,2);
    end
    x_ax_string=nodes(sel_nodes);
elseif strcmp(var_type,'state')
    scan_values=stat_sol_states_lhs_parscan;
    if isempty(sel_nodes)
        sel_nodes=1:size(scan_values,2);
    end
    x_ax_string=arrayfun(@(x) strcat('state #',num2str(x)), sel_nodes,'un',0);
else
    error('<var_type> must be ''state'' or ''node''.')
end

if isempty(sobol_sensit_index)

if ~isempty(sample_size)
    M=sample_size;
else
    M=size(all_par_vals_lhs_subset,1)/2; 
end

A=all_par_vals_lhs_subset(1:M,:); B=all_par_vals_lhs_subset(M+1:2*M,:);
tau_i=zeros(size(all_par_vals_lhs_subset,2),numel(sel_nodes));
% transition_rates_table=ones(2,numel(nodes));

if ~isempty(disp_freq)
  disp_var=disp_freq;
else
  disp_var=[];
end

    for k=1:size(all_par_vals_lhs_subset,2)
    %%%%
    fprintf(strcat('\n','recalculating variance for parameter #',num2str(k),'\n\n'))
    B_i = A; B_i(:,k) = B(:,k);
        % rerun calculations
      if strcmp(var_type,'node')
          [f_B_i,~]=fcn_calc_paramsample_table(B_i,scan_params_sobol,scan_params_up_down_sobol,transition_rates_table,stg_cell,x0,disp_var);
      elseif strcmp(var_type,'state')
          [~,f_B_i]=fcn_calc_paramsample_table(B_i,scan_params_sobol,scan_params_up_down_sobol,transition_rates_table,stg_cell,x0,disp_var);
      end
      
        % sensit index
        for varcount=1:numel(sel_nodes)
            diff_fA_fBi=( scan_values(1:M,sel_nodes(varcount)) - f_B_i(:,sel_nodes(varcount)) );
            if sum(isnan(diff_fA_fBi))>0
                disp(strcat(num2str(100*sum(isnan(diff_fA_fBi))/numel(diff_fA_fBi)),'% nans in solutions!!'));
                diff_fA_fBi=diff_fA_fBi(~isnan(diff_fA_fBi));
            end
            tau_i(k,varcount)=(diff_fA_fBi'*diff_fA_fBi)/( 2*M*nanvar(scan_values(1:M,sel_nodes(varcount))) );
        end
    
    end

else
    tau_i=sobol_sensit_index;
end

if ~isempty(plot_settings)
num_size_plot=plot_settings(1); 
if length(plot_settings)>=5 && ~isnan(plot_settings(4)) && ~isnan(plot_settings(5))
    min_col_val=plot_settings(4); maxval_color=plot_settings(5); % 1.05*max(tau_i(:))
else
    min_col_val=0; maxval_color=1; % 1.05*max(tau_i(:))
end

heatmap(tau_i,x_ax_string,predictor_names,'%0.2f','TickAngle',plot_settings(6),'Colormap','redblue',...
    'MinColorValue',min_col_val,'MaxColorValue',maxval_color,'GridLines','-',...
    'FontSize',num_size_plot,'ShowAllTicks',true,'colorbar',true,'NaNColor',[0 0 0])
set(gca,'FontSize',plot_settings(2)); title('Sobol total sensitivity index', 'Fontweight','normal', 'FontSize',plot_settings(3));
end