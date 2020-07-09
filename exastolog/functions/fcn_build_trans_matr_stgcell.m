function [A_sparse_fast,K_sparse]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,kin_matr_flag)

dim_matr=2^size(transition_rates_table,2); n_nodes=size(transition_rates_table,2);
% state_transitions_inds=[trans_source_states_mat, trans_target_states_mat, cell2mat(node_inds), up_down_inds_arr];
    
% trans_source_states_mat=stg_table(:,1); trans_target_states_mat=stg_table(:,2); up_down_inds_arr=stg_table(:,4); 
rate_inds_numel=cell2mat(cellfun(@(x) numel(x),stg_cell,'un',0));
rate_inds_matr=cell2mat(arrayfun(@(n) cell2mat(vertcat(arrayfun(@(x) repelem([x n],rate_inds_numel(n,x),1), 1:size(rate_inds_numel,2),'un',0))'),...
    1:size(rate_inds_numel,1),'un',0)');
rate_inds=(rate_inds_matr(:,1)-1)*2 + rate_inds_matr(:,2);

% up_trans=vertcat(stg_cell{1,:}); down_trans=vertcat(stg_cell{2,:});

B=sparse([vertcat(stg_cell{1,:}); vertcat(stg_cell{2,:})], ...
         [vertcat(stg_cell{1,:}); vertcat(stg_cell{2,:})] + [repelem(1,numel(vertcat(stg_cell{1,:})),1);...
         repelem(-1,numel(vertcat(stg_cell{2,:})),1)].*(2.^(n_nodes - rate_inds_matr(:,1))), ...
    transition_rates_table(rate_inds)/sum(transition_rates_table(:)),dim_matr,dim_matr);
% clearvars stg_cell

A_sparse_fast = B + (speye(size(B)) - diag(sum(B,2)));

if ~isempty(kin_matr_flag)
    K_sparse=(transpose(A_sparse_fast) - speye(size(A_sparse_fast)) )*sum(transition_rates_table(:));
else
    K_sparse=[];
end
