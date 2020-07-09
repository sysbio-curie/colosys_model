function []=fcn_write_logicrules(nodes,rules,truth_table_filename)

rules_rewritten=cell(1,numel(nodes));

all_rules_strings=cellfun(@(x) strsplit(x,{'|','&','!','~','(',')'}), strrep(rules,' ',''), 'UniformOutput',false);
all_rules_strings=horzcat(all_rules_strings{:}); nodes_in_rules=unique(strtrim(all_rules_strings(cellfun(@(x) ~isempty(x), all_rules_strings))));
if sum(ismember(nodes_in_rules,nodes))~=numel(nodes_in_rules)
    error('unknown node names in rules, check the logical rules again')
end

for k=1:numel(rules)

% replace exclamation marks (negation in BoolNet) with '~' (negation in MATLAB)
if ~isempty(strfind(rules{k},'!'))
    rules{k}=strrep(rules{k},'!','~');
end

node_match_cell=regexp(rules{k},'\w*','match'); 
% check if 'rules' have the node names of 'nodes'
% if sum(ismember(node_match_cell,nodes))~=numel(node_match_cell)
%     disp( strcat( node_match_cell(~ismember(node_match_cell,nodes)), '' ) )
% end

% [tf,loc]=ismember(nodes,node_match_cell); [~,p]=sort(loc(tf)); idx=find(tf); idx=idx(p);
[~,b]=ismember(node_match_cell,nodes); inds_print=arrayfun(@num2str, b, 'UniformOutput', false);
node_repl_cell=strcat('list_binary_states(:,',inds_print,')');
rules_rewritten{k}=regexprep(rules{k},node_match_cell,node_repl_cell);

if isempty(rules{k})
    rules_rewritten{k}='zeros(2^n,1)'; % zeros(2^n,1)
end

end

% print to a file
if exist(truth_table_filename,'file')==2
    disp('file exists - overwritten')
    delete(truth_table_filename)
end
fileID1 = fopen(truth_table_filename,'w'); 
fprintf(fileID1,strcat('function update_matrix = ',strrep(truth_table_filename,'.m',''),'(nodes)','\n\n', ...
    'n=numel(nodes); list_binary_states=fliplr(rem(floor((0:((2^n)-1)).''*pow2(0:-1:-n+1)),2));','\n\n', ...
    'update_matrix = ['));

for k=1:numel(rules_rewritten)
    if k<numel(rules_rewritten)
        line_end=',... '; formatSpec_string=strcat('%s','\n');
    else
        line_end=''; formatSpec_string=strcat('%s');
    end
    % disp(strcat(rules_rewritten{k},line_end))
    fprintf(fileID1,formatSpec_string,strcat(rules_rewritten{k},line_end));
end
fprintf(fileID1,'];'); fclose(fileID1);
