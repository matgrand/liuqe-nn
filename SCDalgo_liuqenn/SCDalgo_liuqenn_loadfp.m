function FP = SCDalgo_liuqenn_loadfp(obj)

% (here?) compile the mex function -> will generate net_forward2_mex.mex* file
% build_net_forward2_mex;

FP = struct();

% sampling time
FP.dt = 1e-4;

% input dimensions
FP.nBm = int32(38);  % magnetic probes
FP.nFf = int32(38);  % flux loops
FP.nFt = int32(1);    
FP.nIa = int32(19);  % coils
FP.nIp = int32(1);   % Ip
FP.nIu = int32(38);  
FP.nrBt = int32(1);   

FP.nphys = 96;

% FP.nphys = FP.nBm + FP.nFf + FP.nFt + FP.nIa + FP.nIp + FP.nIu + FP.nrBt;
% assert (FP.nphys == 136);

% number of points (could be a tunable par? Actually it's max is fixed at compilation time)
FP.npts = int32(25);

end