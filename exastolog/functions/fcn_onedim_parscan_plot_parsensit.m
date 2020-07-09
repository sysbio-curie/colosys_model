function [resp_coeff,scan_pars_sensit,scan_params_sensit_up_down,fig_filename]=fcn_onedim_parscan_plot_parsensit(plot_types,plot_type_options,...
                                         stationary_node_vals_onedimscan,stationary_state_vals_onedimscan,nonzero_states_inds,parscan_matrix,...
                                         nodes,scan_params,scan_params_up_down,sensit_cutoff,param_settings)
% plot_type_flag: heatmap or lineplot
% var_type_flag: states or nodes
% readout_type_flag: variable value or response coefficient
% var_types={'nodes','states'}

z=arrayfun(@(x) plot_types{x}(plot_type_options(x)), 1:numel(plot_type_options),'un',0); z=horzcat(z{:}); 
str_plot_type=cell2mat(strcat(z,'_'));
fig_filename=strcat('onedim_parscan_',str_plot_type,param_settings{4},'_by_vars');
fig_filename
plot_type_flag=z{1}; var_type_flag=z{2};readout_type_flag=z{3};

if strfind(var_type_flag,'node')
    scan_variable=stationary_node_vals_onedimscan; 
elseif strfind(var_type_flag,'state')
    scan_variable=stationary_state_vals_onedimscan;
else
    error('<var_type_flag> should be ''states'' or ''nodes''!');
end

[~,scan_par_inds,~]=fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
fontsize_axes=param_settings{1}; fontsize_title=param_settings{2};
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);

trans_rates_names={strcat('u_',nodes);strcat('d_',nodes)}'; 
trans_rates_names=vertcat(trans_rates_names{:}); trans_rates_names=horzcat(trans_rates_names(:))';
nonzero_states=truth_table_inputs(nonzero_states_inds,:); 

% disp(trans_rates_names)

% ONE STATE on one subplot as fcn of all sensitive params
if strcmp(var_type_flag,'state') || strcmp(var_type_flag,'states')
    % scan_variable = stationary_state_vals_onedimscan;
    if size(scan_variable,3)~=size(nonzero_states,1)
        error('input variable does not have correct dimension (should be ''stationary_state_vals_onedimscan'')')
    end
elseif strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
    % scan_variable = stationary_node_vals_onedimscan;
     if size(scan_variable,3)~=numel(nodes)
        error('input variable does not have correct dimension (should be ''stationary_node_vals_onedimscan'')')
    end
end

% calculate resp coeffs
num_diff_vals=diff(scan_variable,1,2); num_diff_parmatr = diff(parscan_matrix)';
% calc resp coeffs
resp_coeff=zeros(size(scan_variable)); resp_coeff=resp_coeff(:,2:end,:);
% resp_coeff_min_max = [floor(min(resp_coeff(:))*10)/10 ceil(max(resp_coeff(:))*10)/10];
for k=1:size(scan_variable,3)
    p_div_x = parscan_matrix'./scan_variable(:,:,k);
    resp_coeff(:,:,k) = (num_diff_vals(:,:,k)./num_diff_parmatr).*p_div_x(:,2:end);
end

% identify parameters that have an effect on stat vars

if strcmp(readout_type_flag,'sensitivity') || strcmp(readout_type_flag,'sensit')
   sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); 
   sensit_params_table=vertcat(sensit_params_table{:});
elseif strcmp(readout_type_flag,'value') || strcmp(readout_type_flag,'values')
    % if we are plotting values, it is the minimal variation in the value of selected variables
   sensit_params_table=arrayfun(@(x) (max(scan_variable(:,:,x),[],2)-min(scan_variable(:,:,x),[],2))'>sensit_cutoff,1:size(scan_variable,3),'un',0);
   sensit_params_table=vertcat(sensit_params_table{:});
else
    error('sensitivity or values')
end    
sensit_pars=find(sum(sensit_params_table)>0); sensit_vars=find(sum(sensit_params_table,2)>0)'; 
% colors for all parameters
all_colors=distinguishable_colors(numel(sensit_pars));

% take subspace of respcoeffs and state_vars where there is an effect
resp_coeff_sensit_parts=resp_coeff(sensit_pars,:,sensit_vars);
scan_variable_sensit_parts=scan_variable(sum(sensit_params_table)>0,:,sum(sensit_params_table,2)>0);
% sensit_vars

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS

% PLOTs showing the stationary value of variables
if strcmp(readout_type_flag,'var_value') || strcmp(readout_type_flag,'value') || strcmp(readout_type_flag,'values')

nrow=round(sqrt(size(scan_variable_sensit_parts,3))); ncol=nrow;
if nrow*ncol<size(scan_variable_sensit_parts,3)
    nrow=nrow+1;
end
% ONE NODE on one subplot a.a.f of all params

%%% LINEPLOT, VARIABLE VALUE 
if strcmp(plot_type_flag,'lineplot') || strcmp(plot_type_flag,'line')
    
% nrow=round(sqrt(size(scan_variable,3))); ncol=nrow;
% if nrow*ncol<size(scan_variable,3)
%     nrow=nrow+1;
% end

if ~isempty(param_settings{3})
    t_pars=param_settings{3}; [ha,~]=tight_subplot(nrow,ncol,t_pars{1},t_pars{2},t_pars{3});
end

for k=1:size(scan_variable_sensit_parts,3)
    if isempty(param_settings{3})
        subplot(nrow,ncol,k); 
    else
        axes(ha(k));
    end
    
    % subsetting only those params where there is variation
    sensit_pars_indiv_var=find(max(scan_variable_sensit_parts(:,:,k),[],2) - min(scan_variable_sensit_parts(:,:,k),[],2)>sensit_cutoff)';
    % color_vals = all_colors(ismember(sensit_vars_all,sensit_vars_indiv_param),:);
    
    for j=1:numel(sensit_pars_indiv_var)
        
    semilogx(parscan_matrix(:,sensit_pars(sensit_pars_indiv_var(j))), ...
        scan_variable_sensit_parts(sensit_pars_indiv_var(j),:,k)', 'LineWidth',2, 'Color', all_colors(sensit_pars_indiv_var(j),:) ); 
    hold on;
    end
    
 ylim([0 1]); legend(trans_rates_names(scan_par_inds(sensit_pars(sensit_pars_indiv_var))),'Interpreter','none');
 sensit_pars(sensit_pars_indiv_var)
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(strrep(nodes(sensit_vars(k)),'_','\_'), 'FontWeight','normal','FontSize',fontsize_title); ylim([0 1]); 
    else
       title(strcat(num2str(nonzero_states_inds(sensit_vars(k))),' (',strrep(num2str(nonzero_states(sensit_vars(k),:)),' ',''),')'), ...
           'FontWeight','normal','FontSize',fontsize_title); 
    end
    if rem(k,ncol)==1
        ylabel('stat. probability','Fontsize',fontsize_axes)
    end
    
% if k==size(scan_variable_sensit_parts,3)
%     legend(strrep(trans_rates_names(scan_par_inds(sensit_pars)),'_','\_'), 'Location', 'eastoutside');
% end

end

% h_supt=suptitle('stationary value of variables'); set(h_supt,'Fontsize',1.5*fontsize_axes)

%%%%%%%%%%%%%%%% 
%%% HEATMAP, VARIABLE VALUE
elseif strcmp(plot_type_flag,'heatmap') || strcmp(plot_type_flag,'hmap')

if ~isempty(param_settings{3})
    t_pars=param_settings{3};
    [ha,~]=tight_subplot(nrow,ncol,t_pars{1},t_pars{2},t_pars{3});
end

% resp_coeff_sensit_parts = reshape( resp_coeff_sensit_parts, size(resp_coeff_sensit_parts,2)*size(resp_coeff_sensit_parts,3), size(resp_coeff_sensit_parts,1) );
% resp_coeff_sensit_parts: [params, param values, variables]
for k=1:size(scan_variable_sensit_parts,3)
    if k>=size(scan_variable_sensit_parts,3)-ncol+1
        xlabs = trans_rates_names(scan_par_inds(sum(sensit_params_table)>0));
    else
        xlabs =[];
    end
    
    
    if rem(k,ncol)==1
        % ylabel('stationary value', 'FontSize', fontsize_axes); 
        round_vals=round(log10(parscan_matrix(:,sensit_pars(1))),1)';
        yticks_str=arrayfun(@(x) strcat('10e',num2str(round_vals(x))), 1:numel(round_vals), 'un', 0);
    else
        yticks_str=[];
    end

    if isempty(param_settings{3})
        subplot(nrow,ncol,k); 
    else
        axes(ha(k));
    end

var_to_plot = scan_variable_sensit_parts(:,:,k)'; % resp_coeff_sensit_parts
heatmap(var_to_plot,xlabs,yticks_str,[],'TickAngle',90,'Colormap','redblue',...
    'MinColorValue',0,'MaxColorValue',max(var_to_plot(:)),'GridLines','none','FontSize',11,'ShowAllTicks',true);    

graph_props=get(gca); graph_props.XAxis.FontSize=fontsize_axes;
    % max(abs(resp_coeff_min_max))
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title( nodes(sensit_vars(k)),'Interpreter','none', 'FontWeight','normal','FontSize',fontsize_title)
    else
        title(strcat(num2str(nonzero_states_inds(sensit_vars(k))),', [',strrep(num2str(nonzero_states(sensit_vars(k),:)),' ',''),']'), ...
            'FontWeight','normal','FontSize',fontsize_title);
    end

if isempty(param_settings{numel(param_settings)})
    if k==size(scan_variable_sensit_parts,3)
        colorbar('EastOutside');
    end
else
    colorbar('East');
end

end % end of for loop

% h_supt=suptitle('stationary value of variables'); set(h_supt,'Fontsize',1.5*fontsize_axes)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT RESPONSE COEFFICIENTS (sensitivity metric)

elseif strcmp(readout_type_flag,'sensitivity') || strcmp(readout_type_flag,'respcoeff') || strcmp(readout_type_flag,'resp_coeff')

% nrow=ceil(sqrt(size(resp_coeff_sensit_parts,3))); ncol=ceil(sqrt(size(resp_coeff_sensit_parts,3)));
nrow=round(sqrt(size(resp_coeff_sensit_parts,3))); ncol=nrow;
if nrow*ncol<size(resp_coeff_sensit_parts,3)
    nrow=nrow+1;
end

%%%%%%%%%%%%    
% RESP COEFF, LINEPLOT
if strcmp(plot_type_flag,'lineplot') || strcmp(plot_type_flag,'line')
% LINEPLOT w resp coeffs

if ~isempty(param_settings{3})
    t_pars=param_settings{3};
    [ha,~]=tight_subplot(nrow,ncol,t_pars{1},t_pars{2},t_pars{3});
end

for k=1:size(resp_coeff_sensit_parts,3)
%     p_div_x = parscan_matrix'./stationary_node_vals_onedimscan(:,:,k);
%     resp_coeff(:,:,k) = (num_diff_vals(:,:,k)./num_diff_parmatr).*p_div_x(:,2:end);
% PLOT
    if isempty(param_settings{3})
        subplot(nrow,ncol,k); 
    else
        axes(ha(k));
    end

resp_coeff_var=resp_coeff_sensit_parts(:,:,k)';
% lineplot
semilogx(parscan_matrix(2:end,sensit_pars), resp_coeff_var, 'LineWidth',2); 
ylim([floor(min(resp_coeff_sensit_parts(:))*10)/10 ceil(max(resp_coeff_sensit_parts(:))*10)/10])

if k==size(resp_coeff_sensit_parts,3) 
    legend(strrep(trans_rates_names(scan_par_inds(sensit_pars)),'_','\_'), 'Location', 'EastOutside');
end
if rem(k,ncol)==1
    ylabel('resp. coeff.', 'FontSize', fontsize_axes)
end

if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title(nodes(sensit_vars(k)),'Interpreter','none', 'FontWeight','normal','FontSize',fontsize_title)
else
        title(strcat(num2str(nonzero_states_inds(sensit_vars(k))),', ',strrep(num2str(nonzero_states(sensit_vars(k),:)),' ',''),']'),...
            'FontWeight','normal','FontSize',fontsize_title); 
end

end

% h_supt=suptitle('response coefficients'); set(h_supt,'Fontsize',1.5*fontsize_axes)

%%%%%%%%%%%%    
% RESP COEFF, HEATMAP
else
% heatmap
% parameters that have abs(resp_coeff)>cutoff
% sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); 
% sensit_params_table=vertcat(sensit_params_table{:});
resp_coeff_sensit_parts=resp_coeff(sum(sensit_params_table)>0,:,sum(sensit_params_table,2)>0);

sensit_vars=find(sum(sensit_params_table,2)>0)';
% resp_coeff_sensit_parts = reshape( resp_coeff_sensit_parts, size(resp_coeff_sensit_parts,2)*size(resp_coeff_sensit_parts,3), size(resp_coeff_sensit_parts,1) );
% resp_coeff_sensit_parts: [params, param values, variables]
abs_min_max=abs([floor(min(resp_coeff_sensit_parts(:))*10)/10 ceil(max(resp_coeff_sensit_parts(:))*10)/10]);
color_min_max=[-max(abs_min_max) max(abs_min_max)];

if ~isempty(param_settings{3})
    t_pars=param_settings{3};
    [ha,~]=tight_subplot(nrow,ncol,t_pars{1},t_pars{2},t_pars{3});
end


for k=1:size(resp_coeff_sensit_parts,3)
    if k>=size(resp_coeff_sensit_parts,3)-ncol+1
        xlabs = trans_rates_names(scan_par_inds(sum(sensit_params_table)>0));
    else
        xlabs =[];
    end

if isempty(param_settings{3})
        subplot(nrow,ncol,k); 
    else
        axes(ha(k));
end

var_to_plot = resp_coeff_sensit_parts(:,:,k)'; % resp_coeff_sensit_parts

    if rem(k,ncol)==1
        % ylabel('response coeff.','FontSize',fontsize_axes)
        round_vals=round(log10(parscan_matrix(:,sensit_pars(1))),1)';
        yticks_str=arrayfun(@(x) strcat('10e',num2str(round_vals(x))), 1:numel(round_vals), 'un', 0);
    else
        yticks_str=[];
    end
    
heatmap(var_to_plot,xlabs,yticks_str,[],'TickAngle',90,'Colormap','redblue',...
    'MinColorValue',color_min_max(1),'MaxColorValue',color_min_max(2),'GridLines','none','FontSize',11,'ShowAllTicks',true);    
    % max(abs(resp_coeff_min_max))
    
    if strcmp(var_type_flag,'node') || strcmp(var_type_flag,'nodes')
        title( nodes(sensit_vars(k)),'Interpreter','none', 'FontWeight','normal','FontSize',fontsize_title)
    else
        title(strcat(num2str(nonzero_states_inds(sensit_vars(k))),', [',strrep(num2str(nonzero_states(sensit_vars(k),:)),' ',''),']'), ...
            'FontWeight','normal','FontSize',fontsize_title); 
    end

    if k==size(resp_coeff_sensit_parts,3)
        colorbar('EastOutside');
    end

end

% h_supt=suptitle('response coefficients'); set(h_supt,'Fontsize',1.5*fontsize_axes)

end % end of IF gate for response coefficient plots

else % readout type (value or sensitivity)
    error('readout_type_flag should be ''values/var_value'' or ''respcoeff/sensitivity'' ')
end

% sensit_params_table=arrayfun(@(x) max(abs(resp_coeff(:,:,x)'))>sensit_cutoff, 1:size(resp_coeff,3),'un',0); sensit_params_table=vertcat(sensit_params_table{:});
% sensit_pars=find(sum(sensit_params_table)>0);

[scan_par_table,~,~]= fcn_get_trans_rates_tbl_inds(scan_params,scan_params_up_down,nodes);
% sensitive parameters
scan_pars_sensit=unique(scan_par_table(sum(sensit_params_table)>0,1))';
scan_params_sensit_up_down=arrayfun(@(x) scan_par_table(sensit_pars(scan_par_table(sensit_pars,1)==x), 2), scan_pars_sensit, 'un', 0);

if any(cell2mat(arrayfun(@(x) size(scan_params_sensit_up_down{x},1),1:numel(scan_params_sensit_up_down),'un',0))>1)
    scan_params_sensit_up_down=arrayfun(@(x) scan_params_sensit_up_down{x}',1:numel(scan_params_sensit_up_down),'un',0);
end

hold off