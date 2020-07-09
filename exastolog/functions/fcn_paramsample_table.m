function paramsample_table = fcn_paramsample_table(scan_vals,n_pars)

n_scanvals=numel(scan_vals); % n_pars=numel(predictor_names);
paramsample_table=scan_vals(fliplr(rem(floor([0:((n_scanvals^n_pars)-1)].'*n_scanvals.^(0:-1:-n_pars+1)),n_scanvals))+1);
