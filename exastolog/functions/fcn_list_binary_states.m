function [list_binary_states,sel_states]=fcn_list_binary_states(n,linear_inds)

list_binary_states=fliplr(rem(floor((0:((2^n)-1)).'*pow2(0:-1:-n+1)),2));

if ~isempty(linear_inds)
    sel_states=list_binary_states(linear_inds,:);
end