function [x,fval,history] = output_fminsearch(x0)

    history=[];
    options = optimset('OutputFcn', @myoutput);
    [x,fval] = fminsearch(@objfun, x0,options);
        
    function stop = myoutput(x,optimValues,state)
        stop = false;
        if isequal(state,'iter')
          history = [history; x optimValues.fval];
        end
    end
    
    function z = objfun(x)
      z = exp(x(1))*(4*x(1)^2+2*x(2)^2+x(1)*x(2)+2*x(2));
    end

end