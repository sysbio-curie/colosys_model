function [nodes,rules]=fcn_bnet_readin(model_oldfilename)

model_newfilename=strrep(model_oldfilename,'.bnet','.dat');
% rules=bnet_table.factors;
prov_filename=strrep(model_oldfilename,'.bnet','1.bnet');
system(['cp ' model_oldfilename ' ' strrep(model_oldfilename,'.bnet','1.bnet')]); 
system(['mv ' prov_filename ' ' model_newfilename]); bnet_table=readtable(model_newfilename,'Delimiter',','); 
nodes=bnet_table.targets'; rules=bnet_table.factors';
