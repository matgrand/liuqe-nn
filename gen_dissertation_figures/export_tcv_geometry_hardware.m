% Export the TCV magnetic-diagnostic and active-coil geometry required by tcv_geometry.py.
% Run this script from the liuqe-ml repository on LAC, then copy tcv_params/tcv_geometry_hardware.mat to the same relative path in the local repository.
clear; clc;

SHOT = 79742;
script_path = mfilename('fullpath'); repo_dir = fileparts(fileparts(script_path));
output_path = fullfile(repo_dir, 'tcv_params', 'tcv_geometry_hardware.mat');

mdsconnect('tcvdata.epfl.ch');
L = liuqe(SHOT, []);
active_winding = any(L.G.Twa ~= 0, 2);
rm = L.G.rm(:); zm = L.G.zm(:); am = L.G.am(:);
rf = L.G.rf(:); zf = L.G.zf(:);
rw = L.G.rw(active_winding); zw = L.G.zw(active_winding);
ww = L.G.ww(active_winding); hw = L.G.hw(active_winding);
Twa = L.G.Twa(active_winding, :); dima = L.G.dima;
schema_version = 1; shot = SHOT;

save(output_path, 'schema_version', 'shot', 'rm', 'zm', 'am', 'rf', 'zf', 'rw', 'zw', 'ww', 'hw', 'Twa', 'dima');
fprintf('Saved %s\n', output_path);
