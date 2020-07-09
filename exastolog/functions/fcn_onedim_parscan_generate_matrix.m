function par_scan_vals=fcn_onedim_parscan_generate_matrix(scan_params,scan_params_up_down,nodes,sampling_type,parscan_min_max,n_steps)

[~,scan_par_inds,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);

if strcmp(sampling_type,'lin') || strcmp(sampling_type,'linear')
    par_range = linspace(parscan_min_max(1),parscan_min_max(2),n_steps);
elseif strcmp(sampling_type,'log') || strcmp(sampling_type,'lognorm')
    par_range = logspace(log10(parscan_min_max(1)),log10(parscan_min_max(2)),n_steps);
else
error('Check sampling type, it has to be ''log/lognorm'' or ''lin/linear''.')    
end

par_scan_vals=repmat(par_range,numel(scan_par_inds),1)';

