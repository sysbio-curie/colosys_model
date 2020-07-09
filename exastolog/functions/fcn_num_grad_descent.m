function [init_error_table,optim_pars_conv,statsol_parscan,error_conv]=fcn_num_grad_descent(state_var_ind,init_error_table,input_cell,data_param_vals,...
                                                                         init_par_vals,incr_resol,incr_resol_init,error_thresh,step_thresh)

% input_cell = {y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names};
y_data=input_cell{1}; x0=input_cell{2}; stg_cell=input_cell{3}; stg_sorting_cell=input_cell{4}; nodes=input_cell{5};
predictor_names=input_cell{6};
if strcmp(state_var_ind,'var')
    [~,fcn_statsol_values]=fcn_handles_fitting(y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);
    init_error=sum((y_data - fcn_statsol_values(init_par_vals) ).^2); % fcn_statsol_sum_sq_dev(init_par_vals);
else
    [fcn_statsol_states_sum_sq_dev,fcn_statsol_values_states]=fcn_handles_fitting('states',y_data,x0,stg_cell,stg_sorting_cell,nodes,predictor_names);
    init_error=fcn_statsol_states_sum_sq_dev(init_par_vals); 
    % boolean state data is column vector, so we transpose
end

if isempty(init_error_table)
init_error_table=zeros(numel(data_param_vals),2);
for up_down_counter=1:2
for k=1:numel(data_param_vals)
    init_par_vals_mod=init_par_vals; 
if up_down_counter==1; init_par_vals_mod(k)=init_par_vals_mod(k)*(1+incr_resol_init); 
else init_par_vals_mod(k)=init_par_vals_mod(k)*(1-incr_resol_init); end
% disp(init_par_vals_mod-init_par_vals)
% error wrt original guess, if positive, error is growing
    if strcmp(state_var_ind,'var')
        init_error_table(k,up_down_counter)=sum((y_data - fcn_statsol_values(init_par_vals_mod) ).^2) - init_error;
    else
        init_error_table(k,up_down_counter)=fcn_statsol_states_sum_sq_dev(init_par_vals_mod) - init_error;
    end
end
end
end % if empty error table

% fractional change in error: init_error_table./init_error
% in which direction is the error decreasing? 1: param value increased, 2: param value decreased
direction_row=arrayfun(@(x) find(init_error_table(x,:)<0),1:size(init_error_table,1),'un',0);

if sum(cellfun(@(x) numel(x),direction_row)==0)>0
[~,min_inds]=min(init_error_table,[],2);
direction_row{cellfun(@(x) numel(x),direction_row)==0}=min_inds(cellfun(@(x) numel(x),direction_row)==0);
end

if sum(cellfun(@(x) numel(x),direction_row)>1)>0
% [~,b]=min(init_error_table(cellfun(@(x) numel(x),direction_row)>1,:));
[~,min_inds]=min(init_error_table,[],2);
direction_row{cellfun(@(x) numel(x),direction_row)>1}=min_inds(cellfun(@(x) numel(x),direction_row)>1);
end

up_down_inds=cell2mat(direction_row);
% incr_resol=0.02;
incr_vector = repmat(1+incr_resol, size(up_down_inds)); incr_vector(up_down_inds==2)=1-incr_resol;
k=0; error_val=init_error; 
step_truth_val=k<step_thresh; if isempty(step_truth_val); step_truth_val=1; end 

disp('starting incrementing trans. rates by initial gradient')

while error_val>init_error*error_thresh && step_truth_val
    k=k+1;
    
    if k==1; optim_pars=init_par_vals; else optim_pars=optim_pars.*incr_vector; end
    optim_pars_conv(k,:)=optim_pars; 
    if strcmp(state_var_ind,'var')
        statsol_parscan(k,:)=fcn_statsol_values(optim_pars);
        error_val=sum((y_data - statsol_parscan(k,:) ).^2); error_conv(k,1)=error_val; 
    else
        state_vals=fcn_statsol_values_states(optim_pars);
        statsol_parscan(k,:)=full(state_vals(state_vals>0));
        error_val=sum((y_data' - statsol_parscan(k,:) ).^2); error_conv(k,1)=error_val; 
    end
    
    
        
        if k>3  && (error_conv(k,1)>=error_conv(k-1,1) && error_conv(k,1)>=error_conv(k-2,1) )
            break
        end  
    disp(k)
    disp(error_val);
end