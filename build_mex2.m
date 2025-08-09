% build_mex_embedded.m

% 1. Define the input argument types for the inputs that remain.
%    We can handle variable-size r and z as before.

phys_arg = coder.typeof(single(0), [1 136], [0 0]); 

upper_bound = 28*65; % Define a max number of points
r_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);
z_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);

% 2. Generate the MEX function.
%    The arguments list is now shorter.
%    The 'model_weights.mat' file must be on the MATLAB path during this step.
codegen net_forward2 -args {phys_arg, r_arg, z_arg} -report -rowmajor 

disp('Embedded-weight MEX file "net_forward2_mex" has been generated.');