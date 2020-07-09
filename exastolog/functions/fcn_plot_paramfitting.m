function []=fcn_plot_paramfitting(state_var_flag,data_init_optim,T_loss,nodes,sel_nodes,vars_show,thres_ind,plot_settings)

if isempty(sel_nodes); sel_nodes=1:size(data_init_optim,2); end
    
if strcmp(state_var_flag,'var') || strcmp(state_var_flag,'vars')
    xlabels_strings=strrep(nodes(sel_nodes),'_','\_');
else
    xlabels_strings=arrayfun(@(x) strcat('state',num2str(x)),sel_nodes,'un',0);
end

label_ticks_fontsize=plot_settings(1); label_fontsize=plot_settings(2);

if size(T_loss,2)>1 
    error_data=T_loss(1:thres_ind,vars_show);
    init_error=T_loss(1,2); best_error=T_loss(end,2);
else
    error_data=T_loss; 
    if isempty(vars_show); vars_show=2;end
end

if isempty(thres_ind); thres_ind=max(size(T_loss)); init_error=T_loss(1); best_error=T_loss(end); end

fig_subpl1=subplot(1,2,1); 
plot(1:thres_ind, error_data,'LineWidth',4); xlim([1 thres_ind]); if init_error/best_error>30; set(gca,'yscale','log'); end
set(gca,'FontSize',label_ticks_fontsize); legend_strs={'temperature','sum of squared error'}; 
if numel(vars_show)==1; ylabel(legend_strs(2),'FontSize',label_fontsize); else; legend(legend_strs(vars_show)); end
grid on; xlabel('number of iterations','FontSize',label_fontsize); title('fitting convergence','FontSize',label_fontsize*1.2)
set(fig_subpl1,'Position',[0.07 0.11 0.42 0.815]); 
% title('Parameter fitting by simulated annealing','FontSize',22)
fig_subpl2=subplot(1,2,2); 
barplot_gca=barh(data_init_optim(:,sel_nodes)'); set(fig_subpl2,'ytick',1:numel(sel_nodes)); 
legend({'data','initial guess','optimized'},'FontSize',label_fontsize,'Box', 'off')
set(gca,'FontSize',label_ticks_fontsize); 
% xticklabels=get(gca,'xtick'); set(fig_subpl2,'xticklabel',xticklabels,'FontSize',label_ticks_fontsize);
xlabel('stationary probabilities','FontSize',label_fontsize); 
set(fig_subpl2,'yticklabel',xlabels_strings,'FontSize',label_fontsize,'Position',[0.56 0.11 0.42 0.815]);
title('model variables (data, pre-, post-fitting)','FontSize',label_fontsize*1.2)
grid on;
