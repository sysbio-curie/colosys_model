function stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0)

% disp('identifying SCCs')
subnetws=conncomp(digraph(A_sparse,'omitselfloops'),'Type','weak');
num_subnets=length(unique(subnetws)); cell_subgraphs=cell(num_subnets,1); scc_submat_cell=cell(num_subnets,1);

% disp('identifying SCCs in subgraphs')
for counter=1:num_subnets
    submatrix_inds=find(subnetws==counter); cell_subgraphs{counter}=submatrix_inds;
    A_sparse_sub=A_sparse(submatrix_inds,submatrix_inds);
    scc_submat_cell{counter}=conncomp(digraph(A_sparse_sub,'omitselfloops'),'Type','strong');
end

nonempty_subgraphs=find(arrayfun(@(x) sum(x0(subnetws==x))>0, 1:numel(scc_submat_cell) ));
sorted_vertices_cell=cell(numel(nonempty_subgraphs),1);
cyclic_sorted_subgraphs_cell=cell(numel(nonempty_subgraphs),1);

counter=0;
% disp('Sorting nonempty subgraphs')
for k=nonempty_subgraphs
    counter=counter+1; A_sparse_sub=A_sparse(subnetws==k,subnetws==k); 
    % if all SCCs single vertex
    if numel( unique(scc_submat_cell{k}) )==size(A_sparse_sub,1)
        sorted_vertices_cell{counter}=toposort(digraph(A_sparse_sub,'omitselfloops'));
    else % there are multi-vertex SCCs (cycles)
   disp('cycles in STG')
    % if entire graph is one connected component, no reordering needed
    if numel(unique(scc_submat_cell{k}))==1
        sorted_vertices_cell{counter}=find(subnetws==k);
    else
        [vert_topol_sort,term_cycles_ind,~,scc_cell,term_cycle_bounds]=fcn_metagraph_scc(A_sparse_sub);
        cycle_lengths=cellfun(@(x) numel(x),scc_cell); [a,~]=histc(cycle_lengths,1:max(cycle_lengths));
        disp(strcat('cycles of length: ', num2str(unique(cycle_lengths)),' (', num2str(a(a>0)),  ' times)' ))
        cyclic_sorted_subgraphs_cell{counter}={vert_topol_sort,term_cycles_ind,term_cycle_bounds};
    end
    end
end

stg_sorting_cell = {subnetws,scc_submat_cell,nonempty_subgraphs,sorted_vertices_cell,cyclic_sorted_subgraphs_cell};