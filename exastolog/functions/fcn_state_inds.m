function inds=fcn_state_inds(yes_no,n_series_exp,n_isl_exp)

% n_series_exp=n_series_exp-1;
% n_isl_exp=n_isl_exp-1;
% yes_no=yes_no-1;

n_sub=n_series_exp-n_isl_exp+1;

inds=sum([1:2^(n_series_exp-1); ...
    reshape(repmat((1:2^((n_series_exp-1)-(n_sub-1)))+(yes_no-1),2^(n_sub-1),1),1,2^(n_series_exp-1))*2^(n_sub-1)]);