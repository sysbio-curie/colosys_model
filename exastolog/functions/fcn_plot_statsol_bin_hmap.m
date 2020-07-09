function binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,term_verts_cell,nodes,sel_nodes,...
                                                    plot_param_settings,tight_subplot_flag,ranking_flag)

% param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
num_size_plot=plot_param_settings(1); fontsize=plot_param_settings(2); hor_gap=plot_param_settings(3); 
bottom_marg=plot_param_settings(4); left_marg=plot_param_settings(5);
n=numel(nodes); truth_table_inputs=fliplr(rem(floor([0:((2^n)-1)].'*pow2(0:-1:-n+1)),2));

% tight_subplot(Nh, Nw, gap, marg_h, marg_w):
% gap: gaps between the axes in normalized units
% marg_h: margins in height in normalized units
% marg_w: margins in width in normalized units

if issparse(stat_sol);  stat_sol=full(stat_sol); end

% take non-empty
term_vertices_input=term_verts_cell(~cellfun(@(x) isempty(x),term_verts_cell));
nonempty_subgraphs=find(~cellfun(@(x) isempty(x),term_verts_cell))';

% extract
if numel(term_vertices_input)==1
    term_vertices_input=term_vertices_input{1};
end
    
% extract those states with probability above threshold
if ~isempty(prob_thresh)
    cntr=0;
for k=1:numel(term_vertices_input)
    if sum(stat_sol(term_vertices_input{k})>prob_thresh)>0
        cntr=cntr+1;
    if iscell(term_vertices_input{k}); curr_elem=cell2mat(term_vertices_input{k}); else; curr_elem=term_vertices_input{k}; end
%     if numel(term_vertices_input{k})>1 || numel(curr_elem)==1
        term_verts_inds_cell_thresh{cntr}=term_vertices_input{k}(stat_sol(curr_elem)>prob_thresh);
%     else
%         term_verts_inds_cell_thresh{k}={term_vertices_input{k}{1}(stat_sol(curr_elem)>prob_thresh)};
%     end
    end
end
else % isempty(prob_thresh)
    term_verts_inds_cell_thresh=term_vertices_input;
end

if ~isempty(tight_subplot_flag)
    n_states=numel(term_verts_inds_cell_thresh); % sum(arrayfun(@(x) numel(term_verts_inds_cell_thresh{x}), 1:numel(term_verts_inds_cell_thresh)));
    [ha,~]=tight_subplot(n_states,1,[hor_gap hor_gap],[bottom_marg 0.01],[left_marg 0.01]);
end

if isempty(sel_nodes); sel_nodes=1:numel(nodes); end

subplot_counter=0; n_cells=numel(term_verts_inds_cell_thresh);
for k=1:numel(term_verts_inds_cell_thresh)
inds=term_verts_inds_cell_thresh{k}; % disp(k);
if iscell(inds); inds=cell2mat(term_verts_inds_cell_thresh{k}); end
    
    % rank by probability
    if ~isempty(ranking_flag); [~,ranking]=sort(stat_sol(inds)); ranking=flipud(ranking); else; ranking=1:numel(inds); end; n_prec=3;
    % multiple elements in cell
%     if numel(term_verts_inds_cell_thresh{k})>1
%         for inner_c=1:numel(term_verts_inds_cell_thresh{k})
%             subplot_counter=subplot_counter+1;
%                 if subplot_counter==n_states; x_ax_leg=nodes(sel_nodes); else x_ax_leg=[]; end; y_ax_leg=round(stat_sol(inds(ranking)),n_prec); 
%                 if ~isempty(tight_subplot_flag); axes(ha(subplot_counter)); else; subplot(n_states,1,subplot_counter); end
%           binary_heatmap=heatmap(truth_table_inputs(inds(ranking(inner_c)),sel_nodes),...
%               x_ax_leg,strcat(num2str(y_ax_leg(inner_c)), ' (#', num2str(nonempty_subgraphs),')'),... % nonempty_subgraphs(k)
%                 '%0.0f','TickAngle',90,...
%                 'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,'GridLines','-','FontSize',num_size_plot,'ShowAllTicks',true); 
%             set(gca,'FontSize',fontsize)
%         end
%     else
        subplot_counter=subplot_counter+1;
        if k==n_cells; x_ax_leg=nodes(sel_nodes); else x_ax_leg=[]; end; y_ax_leg=round(stat_sol(inds(ranking)),n_prec); 
        if ~isempty(tight_subplot_flag); axes(ha(subplot_counter)); else; subplot(n_states,1,subplot_counter); end
        binary_heatmap=heatmap(truth_table_inputs(inds(ranking),sel_nodes),x_ax_leg,...
            strcat(num2str(y_ax_leg(ranking))), ... % , ' (#', num2str(nonempty_subgraphs),')', nonempty_subgraphs(k)
            '%0.0f','TickAngle',90,'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,'GridLines','-','FontSize',num_size_plot,'ShowAllTicks',true); 
        set(gca,'FontSize',fontsize)
%     end
end

hold off