function [stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,x0_flag)

init_node_vals=[]; 
n=log2(size(x0,1)); 

% truth_table_inputs=rem(floor(up_array.'*pow2(neg_array)),2); 

if (numel(stat_sol)==numel(x0)) && (size(stat_sol,1)==size(x0,1)) && 2^n==numel(stat_sol)
    
    up_array=0:((2^n)-1); neg_array=0:-1:-n+1; 
    stat_states_binary=fliplr(rem(floor(up_array(stat_sol>0).'*pow2(neg_array)),2));
    stationary_node_vals=stat_sol(stat_sol>0)'*stat_states_binary; 

if ~isempty(x0) && ~isempty(x0_flag)
    x0_states_binary = fliplr(rem(floor(up_array(x0>0).'*pow2(neg_array)),2));
    init_node_vals=x0(x0>0)'*x0_states_binary;
else
    init_node_vals=[];
end


else
   disp('initial condition or stationary solution don''t match with each other and/or with number of nodes')
end