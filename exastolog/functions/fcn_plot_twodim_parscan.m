function []=fcn_plot_twodim_parscan(stat_sol_paramsample_table,scanvals,multiscan_pars,multiscan_pars_up_down,nodes,sel_nodes,plot_settings)

up_down_str={'u_','d_'}; label_str=strcat(up_down_str(cell2mat(multiscan_pars_up_down)), nodes(multiscan_pars));
axis_str=arrayfun(@(x) strcat('1e',num2str(x)), round(log10(scanvals),2),'un',0); 
if sum(scanvals==0)>0; axis_str{scanvals==0}='0'; end
if isempty(sel_nodes); sel_nodes=1:numel(nodes); end
n_row_subplot=ceil(sqrt(numel(sel_nodes)));
n_col_subplot=n_row_subplot;
if n_col_subplot*(n_row_subplot-1)>numel(sel_nodes); n_col_subplot=n_col_subplot-1; end

for k=sel_nodes
subplot(n_row_subplot,n_col_subplot,find(k==sel_nodes));
    % y axis is 2nd column of paramsample_table
    heatmap(flipud(reshape(stat_sol_paramsample_table(:,k),numel(scanvals),numel(scanvals))),...
            axis_str,fliplr(axis_str),'%0.2f', ...
            'Colormap','redblue','MinColorValue',0,'MaxColorValue',1,'ShowAllTicks',true,'GridLines','-','FontSize',plot_settings(1));
        set(gca,'XTickLabel',axis_str,'FontSize',plot_settings(1));
        title(nodes(k),'Interpreter','none','FontSize',plot_settings(3)); 
   if find(k==sel_nodes)>n_col_subplot*(n_row_subplot-1)||n_row_subplot==1
       xlabel(label_str{1},'Interpreter','none','FontSize',plot_settings(2)); 
   end
   if (rem(find(k==sel_nodes),n_col_subplot)==1||n_col_subplot==1); ylabel(label_str{2},'Interpreter','none','FontSize',plot_settings(2)); end
end
