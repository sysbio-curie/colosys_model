function [y_data,y_init_pred,init_error]=fcn_param_fitting_data_initguess_error(state_var_flag,x0,stg_cell,data_param_vals,init_par_vals,...
                                                stg_sorting_cell,nodes,predictor_names)

transition_rates_table_optim=fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,data_param_vals);

if strcmp(state_var_flag,'var') || strcmp(state_var_flag,'vars')
y_data=fcn_calc_init_stat_nodevals(x0,...
	split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table_optim,''),...
	stg_sorting_cell,transition_rates_table_optim,x0),'x0');
% initial value of model nodes (with the initial parameter guess)
y_init_pred=fcn_calc_init_stat_nodevals(x0,split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,...
		fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,init_par_vals),''),...
		stg_sorting_cell,transition_rates_table_optim,x0),'');
init_error=sum((y_data-y_init_pred).^2);
else
    
y_data=split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,transition_rates_table_optim,''),stg_sorting_cell,transition_rates_table_optim,x0);
% y_data=full(y_data(y_data>0));
% prediction with init guess
y_init_pred=split_calc_inverse(fcn_build_trans_matr_stgcell(stg_cell,...
		fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,init_par_vals),''),...
		stg_sorting_cell,transition_rates_table_optim,x0);
% y_init_pred=full(y_init_pred(y_init_pred>0));
init_error=sum((y_data-y_init_pred).^2);
    
end