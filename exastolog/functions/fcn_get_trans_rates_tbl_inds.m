function [par_ind_table,trans_rate_scan_inds,predictor_names]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes)

par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];

trans_rate_scan_inds=sub2ind([2 numel(nodes)],par_ind_table(:,2),par_ind_table(:,1));

prefix={'u_','d_'}; 
predictor_names=arrayfun(@(x) strcat(prefix(par_ind_table(x,2)), nodes(par_ind_table(x,1))), 1:size(par_ind_table,1),'un',0);

predictor_names=vertcat(predictor_names{:})'; 