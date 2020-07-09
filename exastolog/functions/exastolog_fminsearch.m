function [x,fval,exitflag,output,history] = exastolog_fminsearch(var_type_flag,y_data,init_par_vals,...
                                                x0,stg_cell,stg_sorting_cell,nodes,predictor_names,fminsearch_opts)

[fcn_statsol_sum_sq_dev,~]=fcn_handles_fitting(var_type_flag,y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);
% fminsearch_opts = optimset('Display','iter','PlotFcns',@optimplotfval);

history=[];
fminsearch_options=optimset('OutputFcn',@myoutput,'Display',fminsearch_opts{1},'MaxFunEvals',fminsearch_opts{2},'MaxIter',fminsearch_opts{3});
[x,fval,exitflag,output]=fminsearch(@(pv) fcn_statsol_sum_sq_dev(pv),init_par_vals,fminsearch_options);
% isempty(fminsearch_options.ActiveConstrTol)

    function stop = myoutput(x,optimValues,state)
        stop = false;
        if isequal(state,'iter')
          history = [history; x optimValues.fval];
        end
    end

%     function z = objfun(x)
%       z = exp(x(1))*(4*x(1)^2+2*x(2)^2+x(1)*x(2)+2*x(2));
%     end

end