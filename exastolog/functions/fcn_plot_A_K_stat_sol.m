% function that plots A (trans matrix), K (kinetic matrix) and solutions
function fcn_plot_A_K_stat_sol(A, nodes, sel_nodes, stat_sol, x0, plot_settings,nonzero_flag)

fontsize=plot_settings(1:3); barwidth_states_val=plot_settings(4);

[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol,'x0');
n=numel(nodes); 
truth_table_inputs=fliplr(rem(floor([0:((2^n)-1)].'*pow2(0:-1:-n+1)),2));

if isempty(sel_nodes)
    sel_nodes=1:numel(nodes); 
end

y_position=0.11;

if ~isempty(A)

min_max_col=plot_settings(5:6); min_col=min_max_col(1); max_col=min_max_col(2);
fontsize_hm=fontsize(1); fontsize_stat_sol=fontsize(2);

% transition matrix
fig_subpl1=subplot(1,3,1); set(fig_subpl1,'Position',[0.032 y_position 0.3 0.82]); % [0.032 0.11 0.32 0.82]
if issparse(A)
    spy(A); set(gca,'FontSize',fontsize_hm)
else
    heatmap(A,1:2^length(nodes),1:2^length(nodes),'%0.2f','TickAngle',90,'MinColorValue',min_col,'MaxColorValue',max_col,'Colormap','redblue',... 
        'GridLines', '-', 'FontSize', fontsize_hm, 'ShowAllTicks', true); 
end
% title 
if full(min(A(:)))==0
    title('A (transition matrix)', 'FontSize', fontsize_stat_sol, 'FontWeight','normal');
else
    title('K (kinetic matrix)', 'FontSize', fontsize_stat_sol, 'FontWeight','normal');
end

%%%%%
% 2nd subplot with states
if ~isempty(nonzero_flag); nnz_vals=find(stat_sol>nonzero_flag)'; state_vals=stat_sol(nnz_vals); else state_vals=stat_sol; end
fig_subpl2=subplot(1,3,2); barplot_object=barh(state_vals,'BarWidth',barwidth_states_val,'FaceColor',[0 0.5 0],'EdgeColor',[0 0.5 0]); 

if ~isempty(nonzero_flag)
    % set(fig_subpl2,'ytick',1:numel(nnz_vals)); set(fig_subpl2,'yticklabel',nnz_vals,'FontSize',fontsize(1)); 
    set(gca,'yticklabel','')
    xlim([0 1.5*max(state_vals)])
    set(fig_subpl2,'Position',[0.34 y_position 0.27 0.82])
    binary_states_cell=arrayfun(@(x) sprintf('%d', truth_table_inputs(nnz_vals(x),:)), 1:numel(nnz_vals),'un',0);
    arrayfun(@(k) text(barplot_object.YData(k)+max(state_vals)/50,barplot_object.XData(k),binary_states_cell{k}, 'FontSize',fontsize(3)),...
        1:numel(binary_states_cell), 'un',0)
else
    set(fig_subpl2,'ytick','') 
end
title('states','FontWeight','normal','FontSize', fontsize(2)); xlabel('stationary probability', 'FontSize', fontsize_hm);

%%%%%
% 3rd subplot with nodes
fig_subpl3=subplot(1,3,3); 
fig_subpl3_vals=flipud([init_node_vals(sel_nodes); stationary_node_vals(sel_nodes)]');
bar_subpl3=barh(fig_subpl3_vals, 'grouped','BarWidth',barwidth_states_val*1.5); grid on; set(gca,'FontSize',fontsize_hm);
xlabel('stationary probability', 'FontSize', fontsize_hm);
title('variables','FontSize',fontsize_stat_sol,'FontWeight','normal'); 
set(fig_subpl3,'Position',[0.72 y_position 0.27 0.815]);
set(gca,'ytick',1:numel(sel_nodes)); set(gca,'YtickLabel',strrep(fliplr(nodes(sel_nodes)),'_',' ')); shift_val=0.4; 
ylim([1-shift_val numel(sel_nodes)+shift_val])
set(bar_subpl3(2),'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);
[~,min_ind]=min(max(fig_subpl3_vals,[],2));
leg=legend(bar_subpl3,{'x_0', 'steady state'},'FontSize',fontsize_stat_sol/1.2,'Location','SouthEast'); % 'Location','SouthEast',
if (min_ind-0.5)/numel(sel_nodes)>y_position
    legend_pos=(min_ind-0.5)/numel(sel_nodes);
else
    legend_pos=y_position;
end
set(leg,'Position',[leg.Position(1) legend_pos leg.Position(3) leg.Position(4)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else % without the kinetic/transition matrix

if ~isempty(nonzero_flag)
    nnz_vals=find(stat_sol>nonzero_flag)'; state_vals=stat_sol(nnz_vals);
else
    state_vals=stat_sol;
end

subplot_x_size=0.44; subplot_x_pos=0.035;

fig_subpl1=subplot(1,2,1); 
barplot_object=barh(state_vals,'BarWidth',barwidth_states_val,'FaceColor',[0 0.5 0],'EdgeColor',[0 0.5 0]); xlim([0 1.25*max(state_vals)])
set(fig_subpl1,'Position',[subplot_x_pos y_position subplot_x_size 0.8]); 
if ~isempty(nonzero_flag)
    set(fig_subpl1,'ytick',1:numel(nnz_vals)); 
    binary_states_cell=arrayfun(@(x) sprintf('%d', truth_table_inputs(nnz_vals(x),:)), 1:numel(nnz_vals),'un',0);
    set(fig_subpl1,'yticklabel',nnz_vals,'FontSize',fontsize(1));
    % bars labeled with states
    arrayfun(@(k) text(barplot_object.YData(k)+max(state_vals)/50,barplot_object.XData(k),binary_states_cell{k}), ...
        1:numel(binary_states_cell), 'un',0)
    % num2cell(nnz_vals)
else
    set(fig_subpl1,'ytick','') 
    % ylim([find(state_vals>disp_lim,1,'first')-0.05*numel(state_vals) find(state_vals>disp_lim,1,'last')+0.05*numel(state_vals)]); 
end

title('states','FontWeight','normal','FontSize', fontsize(2)); xlabel('stationary probability', 'FontSize', fontsize(2)); 
% set(gca, 'FontSize', fontsize)
fig_subpl1.XGrid = 'on'; 

fig_subpl2=subplot(1,2,2); 
if ~isempty(init_node_vals)
    plot_vals=flipud([init_node_vals(sel_nodes); stationary_node_vals(sel_nodes)]');
else
    plot_vals=fliplr(stationary_node_vals(sel_nodes))';
end
bar1=barh(plot_vals,'grouped'); set(gca,'ytick',1:numel(sel_nodes)); set(gca, 'YTickLabel',fliplr(strrep(nodes(sel_nodes),'_',' '))); 
set(fig_subpl2,'Position',[subplot_x_pos+subplot_x_size+0.08 y_position subplot_x_size 0.8]); 
n=size(plot_vals,2); set(bar1(n),'FaceColor',[1 0 0]); set(gca,'FontSize',fontsize(2)/1.5); 
xlim([0 1.05]); ylim([0.5 numel(sel_nodes)+0.5]); grid on
title('nodes','FontWeight','normal','FontSize', fontsize(2)); xlabel('stationary probability','FontSize', fontsize(2)); 
if n>1
legend({'x_0', 'steady state'},'Location','SouthEast','FontSize',fontsize(2));
leg=legend(fig_subpl2,{'x_0', 'steady state'},'FontSize',fontsize(2)/1.2,'Location','SouthEast'); % 'Location','SouthEast',
[~,min_ind]=min(max(plot_vals,[],2));
if (min_ind-0.5)/numel(sel_nodes)>y_position
    legend_pos=(min_ind-0.5)/numel(sel_nodes);
else
    legend_pos=y_position;
end

set(leg,'Position',[leg.Position(1) legend_pos leg.Position(3) leg.Position(4)])
end
    
end