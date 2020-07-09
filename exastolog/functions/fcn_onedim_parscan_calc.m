function [stationary_state_vals_onedimscan,stationary_node_vals_scan,stationary_state_inds_scan]=...
    fcn_onedim_parscan_calc(stg_cell,transition_rates_table,x0,nodes,parscan_matrix,scan_params,scan_params_up_down)
                                                                                            
[~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);

if size(parscan_matrix,2)~=numel(scan_par_inds)
	error('parscan_matrix and scan_par_inds don''t have the same dimension, check them again!')
end

A_dim=2^numel(nodes);
% states corresponding to the transition rates
scan_size=size(parscan_matrix,2);

% trans_matr_inds=arrayfun(@(x) ismember(2*(stg_table(:,3)-1)+stg_table(:,4),x),scan_par_inds,'un',0); % 0.8 sec
% [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table,'');
[A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,'');
stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);

[stat_sol,~,~]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0);
stationary_node_vals_scan=zeros(numel(scan_par_inds),size(parscan_matrix,1),numel(nodes)); % transition_rates_table_mod=transition_rates_table;
stationary_state_vals_onedimscan=zeros( numel(scan_par_inds),size(parscan_matrix,1), sum(stat_sol>0) );
% index of nonzero states
stationary_state_inds_scan = cell(numel(scan_par_inds),size(parscan_matrix,1));
% stationary_state_inds_scan{k,val_counter}=find(stat_sol>0);

n_resol = numel(parscan_matrix(:,1));

% can be parallelized by <parfor>
for k=1:scan_size
    disp(strcat(num2str(round(100*k/scan_size)),'% done'))
    % i_row=stg_table(trans_matr_inds{k},1); j_col=stg_table(trans_matr_inds{k},2); seq_inds=i_row+(j_col-1)*A_dim;
    i_row=stg_cell{scan_par_inds(k)}; 
    if rem(scan_par_inds(k),2)>0; j_col=i_row+2^(numel(nodes) - ceil(scan_par_inds(k)/2)); 
    else; j_col=i_row-2^(numel(nodes) - ceil(scan_par_inds(k)/2)); 
    end
    seq_inds=i_row+(j_col-1)*A_dim;
    % [repelem(1,numel(up_trans),1);repelem(-1,numel(down_trans),1)].*(2.^(n_nodes - rate_inds_matr(:,1)))
    
    for val_counter=1:n_resol
        
        transition_rates_table_mod=transition_rates_table; 
        transition_rates_table_mod(scan_par_inds(k))=parscan_matrix(val_counter,k);
        trans_rate_normalized=transition_rates_table_mod(scan_par_inds(k))/sum(transition_rates_table_mod(:));
        norm_factor = sum(transition_rates_table_mod(:))/sum(transition_rates_table(:));
        % [A_sparse,~]=fcn_build_trans_matr(stg_table,transition_rates_table_mod,'');
        A_sparse_mod=A_sparse/norm_factor; 
        % diagonal to 0, it'll be recalculated
        A_sparse_mod(1:(A_dim+1):numel(A_sparse_mod))=0;
        % reassign relevant trans rates
        A_sparse_mod(seq_inds)=trans_rate_normalized; % 0.05-0.1 sec for 20 nodes
        % diagonal has to be recalculated
        A_sparse_mod = A_sparse_mod + speye(size(A_sparse_mod)) - diag(sum(A_sparse_mod,2)); % 0.33sec
        [stat_sol,~,~]=split_calc_inverse(A_sparse_mod,stg_sorting_cell,transition_rates_table_mod,x0);
        % sols per node
        [stationary_node_vals,~]=fcn_calc_init_stat_nodevals(x0,stat_sol,'');
        stationary_node_vals_scan(k,val_counter,:)=stationary_node_vals;
        stationary_state_vals_onedimscan(k,val_counter,:)=stat_sol(stat_sol>0);
        stationary_state_inds_scan{k,val_counter}=find(stat_sol>0);
    end
end