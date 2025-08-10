phys_arg = coder.typeof(single(0), [1 136], [0 0]); 

upper_bound = 28*65; % Define a max number of points
% upper_bound = 25; % Define a max number of points
r_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);
z_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);

codegen net_forward2 -args {phys_arg, r_arg, z_arg} -report -rowmajor 

disp('Embedded-weight MEX file "net_forward2_mex" has been generated.');