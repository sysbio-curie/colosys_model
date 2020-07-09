function stg_cell=fcn_build_stg_cell(truth_table_filename,nodes)

n_nodes=numel(nodes); 

fcn_name_str=strrep(truth_table_filename,'.m',''); fcn_name=str2func(fcn_name_str); 
update_table=sparse(fcn_name(nodes)); 

up_trans_source=arrayfun(@(x) intersect(find(update_table(:,x)),fcn_state_inds(0,n_nodes,x)), 1:size(update_table,2),'un',0);
down_trans_source=arrayfun(@(x) intersect(find(~update_table(:,x)), fcn_state_inds(1,n_nodes,x)), 1:size(update_table,2),'un',0);


stg_cell=[up_trans_source;down_trans_source];

% down_trans_target=arrayfun(@(x) [[down_trans_source{x}]-2^(x-1) repmat([x 2],numel(down_trans_source{x}),1)], 1:numel(down_trans_source), 'un',0);
% up_trans_target=arrayfun(@(x) [[up_trans_source{x}]+2^(x-1) repmat([x 1],numel(up_trans_source{x}),1)], 1:numel(up_trans_source),'un',0);
% 
% state_transitions_inds=[ [cell2mat(vertcat(down_trans_source(:))); cell2mat(vertcat(up_trans_source(:)))] ...
%            [cell2mat(vertcat(down_trans_target(:))); cell2mat(vertcat(up_trans_target(:)))] ];
