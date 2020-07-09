function l0_blocks=fcn_left_kernel(K_sp_sub_reord,r0_blocks,dim_matr)

disp('constructing left kernel')

dim_kernel=sum(~ismember(sum(r0_blocks,2),0)); 
term_block_inds = (dim_matr-dim_kernel+1):dim_matr; nonterm_block_inds = 1:dim_matr-dim_kernel;
colnum_r_null_array=1:size(r0_blocks,2); 
% left kernel
size_r0_blocks=size(r0_blocks);
l0_blocks=transpose(sparse(size_r0_blocks(1),size_r0_blocks(2))); l0_blocks(transpose(~ismember(r0_blocks,0)))=1; 

% if isa(r0_blocks,'sym'); l0_blocks=sym(l0_blocks); end

X_block=-l0_blocks(colnum_r_null_array,term_block_inds)*K_sp_sub_reord(term_block_inds,nonterm_block_inds)/(K_sp_sub_reord(nonterm_block_inds,nonterm_block_inds));
l0_blocks(colnum_r_null_array,nonterm_block_inds)=X_block;