function output = objfcn_perturb_exps(var_type_flag,fit_vars,stg_cell,nodes,perturb_nodes,...
                                distrib_types,distr_typ_sel,dom_prob,...
                                initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,transition_rates_table)

% initial_fixed_nodes={'cc','KRAS','DSB','CDC25B','g2m_trans','cell_death','CHEK1','MAPKAPK2'}; 
% [chk1i,mk2i]={[0,0],[1,0],[0,1],[1,1]}
% inhib_combs=1-[0,0;1,0;0,1;1,1]; initial_fixed_nodes_vals=[1,1,1,0,0,0];
initial_fixed_nodes_vals_matr=[repmat(initial_fixed_nodes_vals,size(inhib_combs,1),1),inhib_combs]; 
stationary_node_vals_perturbs=zeros(size(inhib_combs,1),numel(nodes)); 

for k=1:size(initial_fixed_nodes_vals_matr,1)
% then we generate the table of transition rates: first row is the 'up'rates, second row 'down' rates, in the order of 'nodes'
    transition_rates_table(1,ismember(nodes,perturb_nodes))=inhib_combs(k,:);
    [A_sparse,~]=fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table,'');
    x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals_matr(k,:),dom_prob,nodes,distrib_types{distr_typ_sel},'');
% Topological sorting of STG and identification of cycles
    stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0);
% calculate stationary solution
    [stat_sol,~,~]=split_calc_inverse(A_sparse,stg_sorting_cell,transition_rates_table,x0);
% by model variables
    stationary_state_vals_perturbs(k,:)=full(stat_sol(stat_sol>0))'; 
    stationary_node_vals_perturbs(k,:)=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0'); 
end

if strcmp(var_type_flag,'var') || strcmp(var_type_flag,'vars')
    output=reshape(stationary_node_vals_perturbs(:,fit_vars),1,size(inhib_combs,1)*numel(fit_vars));
else
    output=reshape(stationary_state_vals_perturbs(:,fit_vars),1,size(inhib_combs,1)*numel(fit_vars));
end