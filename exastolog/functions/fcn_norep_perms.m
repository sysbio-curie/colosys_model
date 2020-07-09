function var_combs_unique = fcn_norep_perms(sel_nodes)

% var_combs=permn(sel_nodes,2); % var_combs = unique(z(:,1:2),'rows'); 
% var_combs_unique=zeros(0.5*numel(sel_nodes)*(numel(sel_nodes)-1),2);
% counter=0;
% for k=1:size(var_combs,1)
%     if all(~ismember(var_combs(1:k-1,[2 1]),var_combs(k,:),'rows')) && var_combs(k,1)~=var_combs(k,2)
%         counter=counter+1;
%         var_combs_unique(counter,:)=var_combs(k,:);
%     end
% end
%

perm_cell=cell(numel(sel_nodes)-1,1);

for k=1:numel(sel_nodes)-1
    val_arr=(sel_nodes(k)+1):sel_nodes(end);
    perm_cell{k}=[repmat(sel_nodes(k),numel(val_arr),1) val_arr'];
end

if all(diff(sel_nodes)==1)
    var_combs_unique = vertcat(perm_cell{:});
else
    var_combs_unique = vertcat(perm_cell{:});
    var_combs_unique=var_combs_unique(ismember(var_combs_unique(:,1),sel_nodes) & ismember(var_combs_unique(:,2),sel_nodes),:);
end