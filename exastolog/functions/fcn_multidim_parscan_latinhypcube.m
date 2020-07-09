function [all_par_vals_lhs,stat_sol_lhs_parscan,stat_sol_states_lhs_parscan]=fcn_multidim_parscan_latinhypcube(min_mean_val,max_stdev_val,sampling_type,...
                                            lhs_scan_dim,scan_params,scan_params_up_down,...
                                            transition_rates_table,stg_cell,x0,nodes)

% lhs_scan_dim=3e3; 

% sampling_type='logunif'; % sampling_type='lognorm' (logarithmically spaced normal distrib) or 
% 'linear' (linearly spaced uniform distribution) or 
% 'logunif' (logarithmically spaced uniform distribution)
n_pars=sum(cellfun(@(x) numel(x), scan_params_up_down));

if numel(min_mean_val)==1 && numel(max_stdev_val)
if strcmp(sampling_type,'lognorm')
    par_min_means=repmat(min_mean_val,1,n_pars); max_stdev_vect=repmat(max_stdev_val,1,n_pars); 
elseif strcmp(sampling_type,'linear')
    % if uniform sampling, then the parameters are the minimum and maximum values of the interval in which we uniformly sample
    par_min_means=repmat(min_mean_val,1,n_pars); max_stdev_vect=repmat(max_stdev_val,1,n_pars); % exponents for log-uniform LHS sampling
elseif  strcmp(sampling_type,'logunif')
    par_min_means=repmat(min_mean_val,1,n_pars); max_stdev_vect=repmat(max_stdev_val,1,n_pars); % exponents for log-uniform LHS sampling
else
    error('Choose for ''lognorm'', ''logunif'' or ''linear'' <sampling_type>')
end
else
    if numel(cell2mat(scan_params_up_down(:)))==numel(min_mean_val) && numel(cell2mat(scan_params_up_down(:)))==numel(max_stdev_val)
    par_min_means=min_mean_val; max_stdev_vect=max_stdev_val;
    else
        error('<par_min_means> and <max_stdev_vect> have to have same length as number of parameters (elements in <scan_params_up_down>)')
    end
end
                                  
% latin hypercube sampling: random sampling within the given parameter ranges
% normal distrib 
n_pars=sum(cellfun(@(x) numel(x), scan_params_up_down));
% par_means=ones(1,n_pars); stdev_vect=ones(1,n_pars)*0.1; 
covar_matr=diag(max_stdev_vect); % corr_val=0;
% normally distributed values
if strcmp(sampling_type,'lognormal') || strcmp(sampling_type,'lognorm') % || strcmp(sampling_type,'log')
    all_par_vals_lhs=10.^lhsnorm(par_min_means,covar_matr,lhs_scan_dim);
elseif strcmp(sampling_type,'log') || strcmp(sampling_type,'logunif') 
    par_max=max_stdev_vect; par_min=par_min_means; 
    all_par_vals_lhs=10.^( lhsdesign(lhs_scan_dim,n_pars).*repmat(par_max-par_min,lhs_scan_dim,1) + repmat(par_min,lhs_scan_dim,1));
% uniformly distributed between 0 and 1, in n intervals
elseif strcmp(sampling_type,'linear') || strcmp(sampling_type,'lin') || strcmp(sampling_type,'linunif')
    par_max=max_stdev_vect; par_min=par_min_means; 
    all_par_vals_lhs=lhsdesign(lhs_scan_dim,n_pars).*repmat(par_max,lhs_scan_dim,1) + repmat(par_min,lhs_scan_dim,1);
end

par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];
trans_rate_scan_inds=sub2ind(size(transition_rates_table),par_ind_table(:,2),par_ind_table(:,1));

transition_rates_table_mod=transition_rates_table;
stat_sol_lhs_parscan=zeros(size(all_par_vals_lhs,1),numel(nodes));

% [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table_mod,'');
[A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table_mod,'');
stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);
[stat_sol,~,~]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table_mod,x0);
% nonzero_states_inds=find(stat_sol>0);
stat_sol_states_lhs_parscan=zeros(size(all_par_vals_lhs,1),sum(stat_sol>0));
% stat_sol_states_lhs_parscan_cell=cell(size(all_par_vals_lhs,1),2);

disp(strcat('dimension of parameter scan:',{' '},num2str(size(all_par_vals_lhs,1)),...
    {' '},'parameter sets of', {' '},num2str(size(all_par_vals_lhs,2)),{' '},'parameters.'))

%%%%%%%%%%%%%%
%%%% START LOOP
% cell_counter=0;
% can be parallelized by <parfor>
for k=1:lhs_scan_dim

transition_rates_table_mod = transition_rates_table;
transition_rates_table_mod(trans_rate_scan_inds)=all_par_vals_lhs(k,:);
[A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table_mod,'');
[stat_sol,~,~]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table_mod,x0);
[stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol,'');
stat_sol_lhs_parscan(k,:)=stationary_node_vals;
nonzero_states=stat_sol(stat_sol>0);
stat_sol_states_lhs_parscan(k,:)=nonzero_states;
% if ~isequal(nonzero_states_inds,find(stat_sol>0))
%     cell_counter=cell_counter+1;
%     stat_sol_states_lhs_parscan_cell{k,1}=nonzero_states; 
%     stat_sol_states_lhs_parscan_cell{k,2}=find(stat_sol>0);
% end

if rem(100*k/lhs_scan_dim,1)==0
    disp(strcat(num2str(round(100*k/lhs_scan_dim)),'% done'))
end

end

% if cell_counter==0
%     stat_sol_states_lhs_parscan_cell=nonzero_states_inds';
% end
