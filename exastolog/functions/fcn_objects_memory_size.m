function sel_objects_table=fcn_objects_memory_size(objects_mem,size_limit_mb)

sel_objects_table=cell2table(arrayfun(@(k) {objects_mem(k).name objects_mem(k).bytes/1e6}, ...
    find(arrayfun(@(x) objects_mem(x).bytes/1e6>size_limit_mb,1:numel(objects_mem))),'un',0)');
