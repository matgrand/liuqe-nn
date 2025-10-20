% build_mex2.m

% Define input types with variable dimensions for the points (r, z)
phys_arg = coder.typeof(single(0), [1 136], [0 0]); % old version (working)
% phys_arg = coder.typeof(single(0), [1 96], [0 0]); % new version

upper_bound = 28*65; % Define a max number of points
% upper_bound = 25; % Define a max number of points
r_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);
z_arg = coder.typeof(single(0), [1 upper_bound], [0 1]);

disp('Starting MEX file generation...');

codegen net_forward2 -args {phys_arg, r_arg, z_arg} -report -rowmajor 

disp('Success! Embedded-weight MEX file "net_forward2_mex" has been generated.');
