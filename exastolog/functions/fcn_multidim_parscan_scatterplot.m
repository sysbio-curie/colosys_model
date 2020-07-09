function []=fcn_multidim_parscan_scatterplot(var_ind,all_par_vals_lhs,scan_values,scan_params,scan_params_up_down,...
                                                nodes,sampling_type,param_settings ) % xlim_val_log

n_bins=param_settings(1); linewidth_val=param_settings(2); 
par_ind_table=[repelem(scan_params, cellfun(@(x) numel(x),scan_params_up_down))', horzcat(scan_params_up_down{:})'];
n_pars=size(par_ind_table,1);

round_sqrt=round(sqrt(n_pars));
if round_sqrt^2 >= n_pars
    n_row_plot=round_sqrt; n_col_plot=n_row_plot;
else
    n_row_plot=round_sqrt; n_col_plot=round_sqrt+1;
end

up_down_str = {'u_','d_'}; 
xlim_val_log=[min(all_par_vals_lhs(:)) max(all_par_vals_lhs(:))]; % min and max of x-axis of subplots
logrange=log10(xlim_val_log);
if ~isnan(param_settings(3))
    label_fontsize = param_settings(3);
else
    label_fontsize =14;
end

for k=1:n_pars
    subplot(n_row_plot,n_col_plot,k);
    % scatter: parameter-variable value
    scatter(all_par_vals_lhs(:,k), scan_values(:,var_ind)); 
    if max(scan_values(:,var_ind))>0; ylim([0 max(scan_values(:,var_ind))]); end % xlim(xlim_val_log); 
    % title( strrep(strcat(nodes{var_ind},', rate: ',up_down_str(par_ind_table(k,2)),nodes(par_ind_table(k,1)) ),'_','\_'), 'FontWeight', 'normal') 
    set(gca,'FontSize', label_fontsize/1.5);
    xlabel( strcat(up_down_str(par_ind_table(k,2)),nodes(par_ind_table(k,1))),'Interpreter','none', 'FontSize', label_fontsize ); % grid on
    % plot means
        % n_bins=round(sqrt(size(all_par_vals_lhs,1)));
        
    if strfind(sampling_type,'log')
        set(gca,'xscale','log'); set(gca,'XTick',10.^round(logrange(1):logrange(end)) ); % round(logrange(1):logrange(end)) (logrange(1):logrange(end))
    % set(gca, 'XMinorTick','off') % , 'YMinorGrid','on'
    % set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
        log_xlabels=strcat('1e',arrayfun(@(x) num2str(x), round(logrange(1):logrange(end)),'un',0)); set(gca,'xticklabel', log_xlabels )
        mean_range_vals=logspace(log10(min(all_par_vals_lhs(:,k))),ceil(log10(max(all_par_vals_lhs(:,k)))),n_bins); 
    else
        mean_range_vals=linspace(min(all_par_vals_lhs(:,k)),ceil(max(all_par_vals_lhs(:,k))),n_bins); 
    end

        [a,b]=histc(all_par_vals_lhs(:,k),mean_range_vals); mean_range_vals=mean_range_vals(a>0);
        var_bin_means=cell2mat(arrayfun(@(x) nanmean( scan_values(b==x,var_ind) ), 1:n_bins,'un',0)); var_bin_means=var_bin_means(a>0);
        hold on; 
        plot(mean_range_vals,var_bin_means,'LineWidth',linewidth_val,'Color',[1 0 0]); 
        hold off;
    if rem(k,n_col_plot)==1 && size(scan_values,2)==numel(nodes)
            ylabel(nodes{var_ind},'Interpreter','none','FontSize',label_fontsize);
    end
end

if size(scan_values,2)==numel(nodes)
        h_supt=suptitle( strcat(nodes{var_ind},' stationary value') ); set(h_supt,'FontSize',label_fontsize*1.5,'FontWeight','normal','Interpreter','none')
elseif size(scan_values,2)~=numel(nodes)
    n_nodes=numel(nodes); % truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
    state_names=param_settings(4:end);
    suptitle(strcat('state #',num2str(var_ind)) ); % 
end

