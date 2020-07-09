function [x,resnorm,exitflag,output,history] = exastolog_lsqnonlin_simult_fitting(var_type_flag,y_data,init_par_vals,lbnds,upbnds,...
                                                stg_cell,nodes,perturb_nodes,predictor_names,... % predictor_names,
                                                fit_vars,...
                                                distrib_types,distr_typ_sel,dom_prob,...
                                                initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,lsqnonlin_opts)

% objective fcn that outputs values of requested variables                            
% objfcn_perturb_exps(var_type_flag,fit_vars,stg_cell,nodes,...
%                 distrib_types,distr_typ_sel,dom_prob,initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,transition_rates_table)

fcn_perturb_exps_optim=@(x)objfcn_perturb_exps(var_type_flag,fit_vars,stg_cell,nodes,perturb_nodes,...
                distrib_types,distr_typ_sel,dom_prob,initial_fixed_nodes,initial_fixed_nodes_vals,inhib_combs,...
                fcn_trans_rates_table(nodes,'uniform',[],[],predictor_names,x));
history=[];

% feed in the function, NOT the sum of squares
fcn_diff_statsol_data = @(p,y) y-fcn_perturb_exps_optim(p);
y=y_data;
% MaxFunctionEvaluations, MaxIterations, FunctionTolerance, OptimalityTolerance, StepTolerance
lsqnonlin_options=optimset('OutputFcn',@myoutput,'Display',...
        lsqnonlin_opts{1},'MaxFunEvals',lsqnonlin_opts{2},'MaxIter',lsqnonlin_opts{3},...
        'TolFun',lsqnonlin_opts{4},'TolX',lsqnonlin_opts{5}); 

[x,resnorm,~,exitflag,output]=lsqnonlin(@(p) fcn_diff_statsol_data(p,y),init_par_vals,lbnds,upbnds,lsqnonlin_options);

    function stop = myoutput(x,optimValues,state)
        stop=false;
        if isequal(state,'iter')
          history=[history; x optimValues.resnorm];
        end
    end

end