function [stat_sol_paramsample_table,stat_sol_states_paramsample_table]=fcn_calc_paramsample_table(paramsample_table,...
                                                                multiscan_pars,multiscan_pars_up_down,...
                                                                transition_rates_table,stg_cell,x0,disp_var)
                                                            
par_ind_table=[repelem(multiscan_pars, cellfun(@(x) numel(x),multiscan_pars_up_down))', horzcat(multiscan_pars_up_down{:})'];
trans_rate_scan_inds=(par_ind_table(:,1)-1)*2 + par_ind_table(:,2);
transition_rates_table_mod=transition_rates_table;
transition_rates_table_prev=transition_rates_table;
n_nodes=size(transition_rates_table,2);

A_dim=2^size(transition_rates_table,2);
% states corresponding to the transition rates
plus_minus_inds = par_ind_table(:,2)'; plus_minus_inds(plus_minus_inds==2)=-1;
row_inds=arrayfun(@(x) stg_cell{par_ind_table(x,2),par_ind_table(x,1)},1:size(par_ind_table,1),'un',0);
col_inds=arrayfun(@(x) stg_cell{par_ind_table(x,2),par_ind_table(x,1)}+plus_minus_inds(x)*2^(n_nodes-par_ind_table(x,1)),1:size(par_ind_table,1),'un',0);
seq_inds=arrayfun(@(x) row_inds{x} + (col_inds{x}-1)*A_dim, 1:numel(col_inds),'un',0); 
all_diag_seq_inds=1:(A_dim+1):A_dim^2;
ind_nums=cellfun(@(x) numel(x),col_inds);
rep_transrate_inds=arrayfun(@(x) repmat(x,1,ind_nums(x)), 1:numel(ind_nums),'un',0);

% n_par=numel(trans_rate_scan_inds);
% [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');

stat_sol_states_paramsample_table=cell(size(paramsample_table,1),1); % zeros(size(all_par_vals_lhs,1),sum(stat_sol>0));
stat_sol_paramsample_table=zeros(size(paramsample_table,1),size(transition_rates_table,2));

parscan_dim=size(paramsample_table,1);
rows_with_zero=unique(paramsample_table>0,'rows'); rows_with_zero=rows_with_zero(sum(rows_with_zero,2)<size(paramsample_table,2),:);
stg_sorting_cell_zeros=cell(1,size(rows_with_zero,1));

disp(strcat('dimension of parameter scan:',{' '},num2str(size(paramsample_table,1)),...
    {' '},'parameter sets of', {' '},num2str(size(paramsample_table,2)),{' '},'parameters.'))

% loop
for k=1:parscan_dim

transition_rates_table_mod(trans_rate_scan_inds)=paramsample_table(k,:);
if k==1
    tic;
    [A_sparse_mod,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table_mod,'');
    stg_sorting_cell=fcn_scc_subgraphs(A_sparse_mod,x0);
    [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell,transition_rates_table_mod,x0);
    toc;
else

transition_rates_table_prev(trans_rate_scan_inds)=paramsample_table(k-1,:);

trans_rate_normalized=transition_rates_table_mod(trans_rate_scan_inds)/sum(transition_rates_table_mod(:));
norm_factor=sum(transition_rates_table_mod(:))/sum(transition_rates_table_prev(:));
% A_sparse_mod=A_sparse_mod/norm_factor; 
% diagonal to 0, it'll be recalculated
A_sparse_mod(all_diag_seq_inds)=0; 
A_sparse_mod=A_sparse_mod/norm_factor;

% replace trans rates that change
A_sparse_mod(cell2mat(vertcat(seq_inds(diff(paramsample_table([k k-1],:))~=0))'))=...
    trans_rate_normalized(cell2mat(vertcat(rep_transrate_inds(diff(paramsample_table([k k-1],:))~=0))));

% recalculate diagonal
% A_sparse_mod=A_sparse_mod + speye(size(A_sparse_mod)) - diag(sum(A_sparse_mod,2)); % 0.33sec
A_sparse_mod = A_sparse_mod + sparse(1:A_dim,1:A_dim,1-sum(A_sparse_mod,2),A_dim,A_dim); % 0.33sec

% calculate solution
if sum(paramsample_table(k,:)==0)>0
    [~,row_match]=ismember(paramsample_table(k,:)>0, rows_with_zero);
    which_zero_ind=row_match(1);
    if isempty(stg_sorting_cell_zeros{which_zero_ind}) 
        stg_sorting_cell_zeros{which_zero_ind}=fcn_scc_subgraphs(A_sparse_mod,x0);
    end
    [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell_zeros{which_zero_ind},transition_rates_table_mod,x0);
else
    [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell,transition_rates_table_mod,x0);
end
end

[stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol,'');
stat_sol_paramsample_table(k,:) = stationary_node_vals;
stat_sol_states_paramsample_table{k} = stat_sol(stat_sol>0);

if ~isempty(disp_var)
if rem(k,disp_var)==0
    disp(strcat(num2str(k),'/',num2str(parscan_dim),', ',num2str(round(100*k/parscan_dim)),'% done'))
end
end

end

if numel(unique(cell2mat(arrayfun(@(x) numel(stat_sol_states_paramsample_table{x}), 1:numel(stat_sol_states_paramsample_table),'un',0))))==1
    stat_sol_states_paramsample_table=cell2mat(stat_sol_states_paramsample_table')';
end