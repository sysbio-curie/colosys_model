function x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_type,plot_flag)

% if initial_fixed_nodes is in a different order than nodes need to rearrange!
[~,b]=ismember(nodes,initial_fixed_nodes);
if ~isequal(sort(b(b>0)),b(b>0))
    initial_fixed_nodes=initial_fixed_nodes(b(b>0)); initial_fixed_nodes_vals=initial_fixed_nodes_vals(b(b>0));
end

n_nodes=numel(nodes); truth_table_inputs=fliplr(rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2));
% define initial values
x0=zeros(2^n_nodes,1); 
% defining a dominant initial state (eg. dom_prob=0.8, ie. 80% probability)
initial_on_nodes_inds=any(cell2mat(arrayfun(@(x) strcmp(nodes,initial_fixed_nodes{x} ), 1:numel(initial_fixed_nodes),'un',0)'));
statespace_decim=sum( bsxfun(@times,truth_table_inputs(:,initial_on_nodes_inds),(2.^fliplr(0:sum(initial_on_nodes_inds)-1)) ), 2);
initial_fixed_nodes_vals_decim=sum(initial_fixed_nodes_vals.*(2.^fliplr(0:numel(initial_fixed_nodes_vals)-1)),2);

% inds_condition=ismember(truth_table_inputs(:,initial_on_nodes_inds),initial_fixed_nodes_vals,'rows');
inds_condition=ismember(statespace_decim,initial_fixed_nodes_vals_decim);

if strfind(distrib_type,'unif')
    x0(inds_condition)=repmat(dom_prob/sum(inds_condition),sum(inds_condition),1); 
    x0(~inds_condition)=repmat((1-dom_prob)/(numel(x0)-sum(inds_condition)),numel(x0)-sum(inds_condition),1); 
elseif strfind(distrib_type,'rand')
    x0(inds_condition)=rand(sum(inds_condition),1); x0=dom_prob*x0/sum(x0);
    x0(~inds_condition)=rand(numel(x0)-sum(inds_condition),1); x0(~inds_condition)=(1-dom_prob)*x0(~inds_condition)/sum(x0(~inds_condition));
else
    error('distrib type shoul be ''uniform'' or ''random''')
end

% rounding precision
n_prec=3;
if round(sum(x0),n_prec)==1
    % disp('sum(x0)=1, OK.')
else
     disp('sum(x0)~=1, something wrong!')
end

if ~isempty(plot_flag)
bar(x0); set(gca,'yscale','log'); xlim([1 2^n_nodes]); % ylim([(1-dom_prob)/2^n_nodes 1])
% subplot(2,1,2); x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes,'broad'); 
% bar(x0); xlim([1 2^13]);set(gca,'yscale','log'); ylim([(1-dom_prob)/2^n_nodes 1])
end