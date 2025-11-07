clc; clear all; close all;

% shot = 87899;
% shot = 87895;
shot = 79742;

DECIMATION = 40; % 4 Decimation factor for the time vector
NDOM = 4; % Number of domains for liuqe (must be >=3 for doublet)
% ITERA = 50; % Max Number of iterations for liuqe NOTE: cant change this bc has some weird values for specific shapes (like 888 for doublets)
ITERFRZ = 20; %20; % Number of Picard iterations before which plasma domain is frozen to ease convergence (inactive if itert>0)
ITERH = 20; %20; % Max Number of FE solver iterations
ITERQ = 10; % iterations for the LCFS

t = liuqe_times(shot); % get times
t = t(1:DECIMATION:end); % decimate times

%% run liuqe
[L, LX, LY] = liuqe(shot,t,'ndom',NDOM,'iterh',ITERH,'iterfrz',ITERFRZ,'iterq',ITERQ);

t = LX.t;
nt = length(t);

%% extract quantities
% Ouputs
Fx = single(LY.Fx); % Plasma poloidal flux map | `(rx,zx,t)` | `[Wb]` |
Iy = single(LY.Iy); % Plasma current density map | `(ry,zy,t)` | `[A/m^2]` |
% Inputs
Bm = single(LX.Bm); % magnetic probe measurements | (*,t) | [T]
Ff = single(LX.Ff); % flux loop poloidal flux | (*,t) | [Wb] 
Ia = single(LX.Ia); % poloidal field coil currents | (*,t) | [A]
Ip = single(LX.Ip);   % plasma current | (t) | [A]
% LCFS
rq = reshape(LY.rq(:,end,:,:), [], size(LY.rq,4)); % (np*nd, nt)
zq = reshape(LY.zq(:,end,:,:), [], size(LY.zq,4)); % (np*nd, nt)


% liuqeplot({{L,LX,LY}},{'MEQ'});

% Quick plot of Fx and Iy at mid time
figure; mid = round(nt/2); sgtitle(sprintf('Shot %d', shot));
mask = (rq(:,mid) ~= 0) & (zq(:,mid) ~= 0);
rq = rq(mask, :);
zq = zq(mask, :);
[rr, zz] = meshgrid(L.rx, L.zx);
subplot(1,2,1); contourf(rr, zz, Fx(:,:,mid), 20, 'LineColor', 'none'); hold on;
scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled'); axis xy; colorbar; title('Fx'); 
[rry, zzy] = meshgrid(L.ry, L.zy);
subplot(1,2,2); contourf(rry, zzy, Iy(:,:,mid), 20, 'LineColor', 'none'); hold on;
scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled'); axis xy; colorbar; title('Iy');


mdsdisconnect; % Disconnect from MDSplus