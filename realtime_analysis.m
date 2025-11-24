clc; clear all; close all;

% load mems from liuqenn
% shot = 88566;
shot = 88567;
% shot = 88564;

% START_TIME = 0.05 % start time for analysis
START_TIME = 0.3 % start time for analysis
STOP_TIME = 1.6 %3.0; % stop time for analysis
dec = 50; % decimation factor


nrt = get_scd_mems(shot,'liuqenn');
Fxn = nrt.liuqenn_Fxn;
Brn = nrt.liuqenn_Brn;
Bzn = nrt.liuqenn_Bzn;
Fx = Fxn.Data; Br = Brn.Data; Bz = Bzn.Data;

figure;
subplot(3,1,1); grid on; plot(Fxn); title('Fxn');
subplot(3,1,2); grid on; plot(Brn); title('Brn');
subplot(3,1,3); grid on; plot(Bzn); title('Bzn');

t = round(Fxn.Time,5); % time vector


% Inputs
fmagin = get_scd_mems(shot,'fmagin');
fmagin = fmagin.fmagin;
% inputs from fmag
Bm = fmagin.Bm.Data;
Bmdot = fmagin.Bmdot.Data;
Ff = fmagin.Ff.Data;
IPHI = fmagin.IPHI.Data;
Ia = fmagin.Ia.Data;
Uf = fmagin.Uf.Data;
Va = fmagin.Va.Data;
enable_exec = fmagin.enable_exec.Data;
feedforwards = fmagin.feedforwards.Data;
references = fmagin.references.Data;
% t = fmagin.t;


rBt = zeros(size(Bm,1),1); % just put zeros for rBt for now

phys = [Bm Ff Ia Uf rBt];
phys_t = round(fmagin.t.Data,5);

%times
valid_t_mask = (t >= START_TIME) & (t <= STOP_TIME);
t = t(valid_t_mask);
Fx = Fx(valid_t_mask, :);
Br = Br(valid_t_mask, :);
Bz = Bz(valid_t_mask, :);
valid_phys_mask = (phys_t >= START_TIME) & (phys_t <= STOP_TIME);
phys = phys(valid_phys_mask, :);
phys_t = phys_t(valid_phys_mask);
assert(max(abs(phys_t - t)) < 1e-8, 'Time vectors do not coincide');



%% analysis


% decimate
Fx = Fx(1:dec:end,:);
Br = Br(1:dec:end,:);
Bz = Bz(1:dec:end,:);
t = t(1:dec:end);
phys = phys(1:dec:end, :);


% figure; grid on; bar(log(abs(phys(1,:)))); title('Bar Plot of phys 1 (RT)');
% fprintf('phys (rt): %s \n', mat2str(phys(1,:),4));
phys_RT = phys;


INTERP_METHOD = 'linear'; % interpolation method for interp2

r = [1.0864 1.0313 0.9357 0.8253 0.7297 0.6746 0.6746 0.7297 0.8253 0.9357 1.0313 1.0864 1.0137 0.9747 0.8805 0.7863 0.7473 0.7863 0.8805 0.9747 0.9182 0.9182 0.8428 0.8428 0.8805];
z = [0.1574  0.4299  0.5873  0.5873  0.4299  0.1574 -0.1574 -0.4299 -0.5873 -0.5873 -0.4299 -0.1574 0. 0.2687 0.38 0.2687 0. -0.2687 -0.38 -0.2687 -0.1075 0.1075 0.1075 -0.1075  0.];

% net inference realtime
Fx_net = Fx; Br_net = Br; Bz_net = Bz;

fprintf('Running LIUQE versions...\n');

% rtliuqe (approx)
[L,LX,LY] = liuqe(shot,t,'itert',1);
[Fx_lrt, Br_lrt, Bz_lrt] = interp_liuqe(r,z, L, LX, LY, INTERP_METHOD);

% lih
[L,LX,LY] = lih('tcv', shot, t, 'debug', 1);
[Fx_lih, Br_lih, Bz_lih] = interp_liuqe(r,z, L, LX, LY, INTERP_METHOD);

% standard liuqe
[L,LX,LY] = liuqe(shot,t);
[Fx_liu, Br_liu, Bz_liu] = interp_liuqe(r,z, L, LX, LY, INTERP_METHOD);

% recalculated net
[Fx_rec, Br_rec, Bz_rec] = liuqe_nn_LX(LX, single(r'), single(z'));


% recalculated on rt phys




% print abs err, std of the err and max err of Fx, Br, Bz between liuqe and liuqenn
e_Fx_nn = abs(Fx_liu - Fx_net);
e_Br_nn = abs(Br_liu - Br_net);
e_Bz_nn = abs(Bz_liu - Bz_net);
fprintf('LIUQENN vs LIUQE:\n');
fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_nn(:)), std(e_Fx_nn(:)), max(e_Fx_nn(:)));
fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_nn(:)), std(e_Br_nn(:)), max(e_Br_nn(:)));
fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_nn(:)), std(e_Bz_nn(:)), max(e_Bz_nn(:)));

% liuqert vs liuqe
e_Fx_lrt = abs(Fx_lrt - Fx_liu);
e_Br_lrt = abs(Br_lrt - Br_liu);
e_Bz_lrt = abs(Bz_lrt - Bz_liu);
fprintf('LIUQERT vs LIUQE:\n');
fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_lrt(:)), std(e_Fx_lrt(:)), max(e_Fx_lrt(:)));
fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_lrt(:)), std(e_Br_lrt(:)), max(e_Br_lrt(:)));
fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_lrt(:)), std(e_Bz_lrt(:)), max(e_Bz_lrt(:)));

% lih vs liuqe
e_Fx_lih = abs(Fx_lih - Fx_liu);
e_Br_lih = abs(Br_lih - Br_liu);
e_Bz_lih = abs(Bz_lih - Bz_liu);
fprintf('LIH vs LIUQE:\n');
fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_lih(:)), std(e_Fx_lih(:)), max(e_Fx_lih(:)));
fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_lih(:)), std(e_Br_lih(:)), max(e_Br_lih(:)));
fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_lih(:)), std(e_Bz_lih(:)), max(e_Bz_lih(:)));

% recalculated net vs liuqe
e_Fx_rec = abs(Fx_rec - Fx_liu);
e_Br_rec = abs(Br_rec - Br_liu);
e_Bz_rec = abs(Bz_rec - Bz_liu);
fprintf('RECALCULATED NET vs LIUQE:\n');
fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_rec(:)), std(e_Fx_rec(:)), max(e_Fx_rec(:)));
fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_rec(:)), std(e_Br_rec(:)), max(e_Br_rec(:)));
fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_rec(:)), std(e_Bz_rec(:)), max(e_Bz_rec(:)));

%plotting
nplt = 5; % number of points to plot

figure;
subplot(3,1,1); grid on; hold on;
plot(t, Fx_net(:,1:nplt:end));
% plot(t, Fx_lrt(:,1:nplt:end), '--');
plot(t, Fx_liu(:,1:nplt:end), ':');
title('Fx'); xlabel('Time (s)');
subplot(3,1,2); grid on; hold on;
plot(t, Br_net(:,1:nplt:end));
% plot(t, Br_lrt(:,1:nplt:end), '--');
plot(t, Br_liu(:,1:nplt:end), ':');
title('Br'); xlabel('Time (s)');
subplot(3,1,3); grid on; hold on;
plot(t, Bz_net(:,1:nplt:end));
% plot(t, Bz_lrt(:,1:nplt:end), '--');
plot(t, Bz_liu(:,1:nplt:end), ':');
title('Bz'); xlabel('Time (s)');

% comparison of the phys inputs
figure;
bm_rt = phys_RT(:,1:38); ff_rt = phys_RT(:,39:76); ia_rt = phys_RT(:,77:95); uf_rt = phys_RT(:,96:133); rbt_rt = phys_RT(:,134);
bm_lx = phys_LX(:,1:38); ff_lx = phys_LX(:,39:76); ia_lx = phys_LX(:,77:95); uf_lx = phys_LX(:,96:133); rbt_lx = phys_LX(:,134);
e_bm=bm_rt-bm_lx; e_ff=ff_rt-ff_lx; e_ia=ia_rt-ia_lx; e_uf=uf_rt-uf_lx; % pure error
e_bm=abs(e_bm); e_ff=abs(e_ff); e_ia=abs(e_ia); e_uf=abs(e_uf); % absolute error

idx = 3;
bmrt=bm_rt(idx,:); ffrt=ff_rt(idx,:); iart=ia_rt(idx,:); ufrt=uf_rt(idx,:); bmlx=bm_lx(idx,:); fflx=ff_lx(idx,:); ialx=ia_lx(idx,:); uflx=uf_lx(idx,:);
% ebm=e_bm(idx,:); eff=e_ff(idx,:); eia=e_ia(idx,:); euf=e_uf(idx,:);

% bmrt=mean(bm_rt,1); ffrt=mean(ff_rt,1); iart=mean(ia_rt,1); ufrt=mean(uf_rt,1); bmlx=mean(bm_lx,1); fflx=mean(ff_lx,1); ialx=mean(ia_lx,1); uflx=mean(uf_lx,1);
ebm=mean(e_bm,1); eff=mean(e_ff,1); eia=mean(e_ia,1); euf=mean(e_uf,1);

nr=4;nc=2;pli=1;
subplot(nr,nc,pli); grid on; hold on; pli=pli+1;
bar(bmrt, 'FaceAlpha',0.5); bar(bmlx, 'FaceAlpha',0.5); title(sprintf('Bm t:%.4f', t(idx))); legend('RT','LX');
subplot(nr,nc,pli); grid on; hold on; pli=pli+1; bar(ebm, 'FaceAlpha',0.5); title('mean(abs(Bm_{rt}-Bm_{lx}))'); 
subplot(nr,nc,pli); grid on; hold on; pli=pli+1;
bar(ffrt, 'FaceAlpha',0.5); bar(fflx, 'FaceAlpha',0.5); title(sprintf('Ff t:%.4f', t(idx))); legend('RT','LX');
subplot(nr,nc,pli); grid on; hold on; pli=pli+1; bar(eff, 'FaceAlpha',0.5); title('mean(abs(Ff_{rt}-Ff_{lx}))'); 
subplot(nr,nc,pli); grid on; hold on; pli=pli+1;
bar(iart, 'FaceAlpha',0.5); bar(ialx, 'FaceAlpha',0.5); title(sprintf('Ia t:%.4f', t(idx))); legend('RT','LX');
subplot(nr,nc,pli); grid on; hold on; pli=pli+1; bar(eia, 'FaceAlpha',0.5); title('mean(abs(Ia_{rt}-Ia_{lx}))');
subplot(nr,nc,pli); grid on; hold on; pli=pli+1;
bar(ufrt, 'FaceAlpha',0.5); bar(uflx, 'FaceAlpha',0.5); title(sprintf('Uf t:%.4f', t(idx))); legend('RT','LX');
subplot(nr,nc,pli); grid on; hold on; bar(euf, 'FaceAlpha',0.5); title('mean(abs(Uf_{rt}-Uf_{lx}))');
% comment on the plot: ia 17 18 19 are probably in wrong order + Uf errors are big
suptitle(sprintf('Shot: %d', shot));



function [fx, br, bz] = interp_liuqe(r, z, L, LX, LY, INTERP_METHOD)
    nt = numel(LY.t);
    np = numel(r);
    fxg = LY.Fx; % Fx on liuqe grid


    i4pirxdzx = 1./(4*pi*L.dzx*L.rx');
    i4pirxdrx = 1./(4*pi*L.drx*L.rx');
    [brg, bzg] = meqBrBz(fxg, i4pirxdzx, i4pirxdrx, L.nzx, L.nrx);

    fx = zeros(nt,np);
    br = zeros(nt,np);
    bz = zeros(nt,np);

    for it = 1:nt
        fx(it,:) = interp2(L.rx, L.zx, fxg(:,:,it), r, z, INTERP_METHOD);
        br(it,:) = interp2(L.rx, L.zx, brg(:,:,it), r, z, INTERP_METHOD);
        bz(it,:) = interp2(L.rx, L.zx, bzg(:,:,it), r, z, INTERP_METHOD);
    end
        

end

% copied from meqpost
function [Br,Bz] = meqBrBz(Fx,i4pirdz,i4pirdr,nz,nr)
    % [Br,Bz] = meqBrBz(Fx,i4pirdz,i4pirdr,nz,nr)
    % Compute Br,Bz fields
    % General version that also accepts time-varying Fx

    [Br,Bz] = deal(zeros(nz,nr,size(Fx,3))); % init
    % Br = -1/(2*pi*R)* dF/dz
    % Central differences dF/dz[i] =  F[i-1] - F[i+1]/(2*dz)
    Br(2:end-1,:,:) = -i4pirdz.* (Fx(3:end,:,:) - Fx(1:end-2,:,:));
    % At grid boundary i, use: dF/dz[i] = (-F(i+2) + 4*F(i+1) - 3*F(i))/(2*dz)
    Br(end,:  ,:) = -i4pirdz          .* (+Fx(end-2,:,:) - 4*Fx(end-1,:,:) + 3*Fx(end,:,:));
    Br(1  ,:  ,:) = -i4pirdz          .* (-Fx(    3,:,:) + 4*Fx(    2,:,:) - 3*Fx(  1,:,:));

    % Bz = 1/(2*pi*R)* dF/dr
    Bz(:,2:end-1,:) =  i4pirdr(2:end-1) .* (+Fx(:,  3:end,:) - Fx(:,1:end-2,:));
    % Same as for Br
    Bz(:,end    ,:) =  i4pirdr(end)     .* (+Fx(:,end-2,:) - 4*Fx(:,end-1,:) + 3*Fx(:,end,:));
    Bz(:,1      ,:) =  i4pirdr(1)       .* (-Fx(:,    3,:) + 4*Fx(:,    2,:) - 3*Fx(:,  1,:));
end
