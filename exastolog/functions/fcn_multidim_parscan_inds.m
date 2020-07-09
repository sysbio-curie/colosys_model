function parscan_inds=fcn_multidim_parscan_inds(n_vals,n_pars)

parscan_inds=fliplr(rem(floor([0:((n_vals^n_pars)-1)].'*n_vals.^(0:-1:-n_pars+1)),n_vals))+1;