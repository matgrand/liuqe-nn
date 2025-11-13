% function to build the MEX file for liuqenn inference

% Define input types with variable dimensions for the points (r, z)
phys_arg = coder.typeof(single(0), [1 134], [0 0]);

% upper_bound = 28*65; % Define a max number of points <- for lih_comparison
upper_bound = 25; % Define a max number of points <- for speed_test
r_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);
z_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);

% LiuqeRT inference mex -> outputs Fx, Br, Bz
codegen liuqenn_inference -args {phys_arg, r_arg, z_arg} -report -rowmajor 

% Iy inference MEX
codegen liuqe_Iy_nn_inference -args {phys_arg, r_arg, z_arg} -report -rowmajor 
