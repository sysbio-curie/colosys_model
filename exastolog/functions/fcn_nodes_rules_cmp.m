function truth_val=fcn_nodes_rules_cmp(nodes,rules)

all_rules_strings=cellfun(@(x) strsplit(x,{'|','&','!','~','(',')'}), strrep(rules,' ',''), 'UniformOutput',false);
all_rules_strings=horzcat(all_rules_strings{:}); 
nodes_in_rules=unique(strtrim(all_rules_strings(cellfun(@(x) ~isempty(x), all_rules_strings))));

if sum(ismember(nodes_in_rules,nodes))==numel(nodes_in_rules)
    truth_val=1;
    disp('Model seems correct: all elements in rules found in nodes list')
else
    truth_val=0;
    bad_node=nodes_in_rules(~ismember(nodes_in_rules,nodes));
    error(strcat('unknown node(s) in rules: {', strjoin(bad_node,', '),'}, check the logical rules again'))
end