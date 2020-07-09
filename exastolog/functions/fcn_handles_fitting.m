function [fcn_statsol_sum_sq_dev,fcn_statsol_values]=fcn_handles_fitting(state_var_flag,y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names)

% if input is 'vars', ie. if we are fitting the stat sol per node
if strcmp(state_var_flag,'var') || strcmp(state_var_flag,'vars')
% create function that calculates sum of squared deviations, composed of different fcns
fcn_statsol_sum_sq_dev=@(x)sum((y_data - fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),''),stg_sorting_cell,...
    fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),x0),'' )).^2);
% function to calculate stationary values
fcn_statsol_values=@(x)fcn_calc_init_stat_nodevals(x0,...
    split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),''),stg_sorting_cell,...
    fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),x0),'' );

else % if input is state prob values

fcn_y_pred_states=@(x)full(split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),''),...
            stg_sorting_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),x0));
% fcn_y_pred_states_nnz=@(x)x(x>0);
fcn_statsol_sum_sq_dev=@(param_vect)sum((y_data - fcn_y_pred_states(param_vect)).^2);
% fcn_statsol_sum_sq_dev=@(param_vect)sum((y_data - fcn_y_pred_states(param_vect)).^2);

% function to calculate stationary values
fcn_statsol_values=@(x)split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),''),...
    stg_sorting_cell,fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x),x0);

end