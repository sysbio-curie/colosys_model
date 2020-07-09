function transition_rates_table=fcn_trans_rates_table(nodes,uniform_or_rand,meanval,sd_val,chosen_rates,chosen_rates_vals)

n=numel(nodes); 

if strcmp(uniform_or_rand,'uniform')
    rate_vals_num=ones(1,2*length(nodes)); % abs(ones(1,2*length(nodes)) + normrnd(0,0.5,1,2*length(nodes)));
elseif strcmp(uniform_or_rand,'random')
    % meanval,sd_val
    rate_vals_num=normrnd(meanval,sd_val,1,2*length(nodes));
    % don't let negative values occur!!
    if any(rate_vals_num<0)
        neg_cnt=0;
        while any(rate_vals_num<0)
            disp('negative value, reassigning')
            rate_vals_num=normrnd(meanval,sd_val,1,2*length(nodes));
            neg_cnt=neg_cnt+1;
            if neg_cnt>100
                break
            end
        end
    end
else
    disp('choose "uniform" or "random" to generate transition rates')
    rate_vals_num=[];
end

if ~isempty(rate_vals_num)

% changing individual transition rate values 
for k=1:numel(chosen_rates)
    split_rate = strsplit(chosen_rates{k},'_'); 
    if numel(split_rate)>2
        node_mod_ind=strjoin(split_rate(2:end),'_');
    else
        node_mod_ind=split_rate{2};
    end
    
    if strcmp(split_rate{1},'d')
        rate_vals_num(find(strcmp(nodes,node_mod_ind))+n)=chosen_rates_vals(k);
    elseif strcmp(split_rate{1},'u')
        % disp(k)
        rate_vals_num(strcmp(nodes,node_mod_ind))=chosen_rates_vals(k); 
    else
        disp('wrong name for transition rate, has to be "u_nodename" or "d_nodename"')
    end
end

% rate_vals_cell=num2cell(rate_vals_num); 
% rate_names (create by strcat)
% [u_cc,u_kras,u_dna_dam,u_chek1,u_mk2,u_atm_atr,u_hr,u_cdc25b,u_g2m_trans,u_cell_death, ...
%     d_cc,d_kras,d_dna_dam,d_chek1,d_mk2,d_atm_atr,d_hr,d_cdc25b,d_g2m_trans,d_cell_death]=deal(rate_vals_cell{:});

transition_rates_table = transpose(reshape(rate_vals_num,length(nodes),2));
end