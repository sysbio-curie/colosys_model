function [x,resnorm,exitflag,output,history] = exastolog_lsqnonlin(var_type_flag,y_data,init_par_vals,lbnds,upbnds,...
                                                x0,stg_cell,stg_sorting_cell,nodes,predictor_names,lsqnonlin_opts)

[~,fcn_statsol_values]=fcn_handles_fitting(var_type_flag,y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);
history=[];
% feed in the function, NOT the sum of squares
fcn_diff_statsol_data = @(p,y) y-fcn_statsol_values(p);
y=y_data;
% MaxFunctionEvaluations, MaxIterations, FunctionTolerance, OptimalityTolerance, StepTolerance
lsqnonlin_options=optimset('OutputFcn',@myoutput,'Display',lsqnonlin_opts{1},'MaxFunEvals',lsqnonlin_opts{2},'MaxIter',lsqnonlin_opts{3},...
    'TolFun',lsqnonlin_opts{4},'TolX',lsqnonlin_opts{5}); % ,'StepTolerance',lsqnonlin_opts{6}

[x,resnorm,~,exitflag,output]=lsqnonlin(@(p) fcn_diff_statsol_data(p,y), init_par_vals,lbnds, upbnds,lsqnonlin_options);
                                   
    function stop = myoutput(x,optimValues,state)
        stop=false;
        if isequal(state,'iter')
          history = [history; x optimValues.resnorm];
        end
    end



end