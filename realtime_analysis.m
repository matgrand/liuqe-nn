clc; clear all; close all;
warning('off', 'all');

% load mems from liuqenn
% shot = 88719; %inverted currents
% shot = 88564; 
shot = 88566;
% shot = 88567;
% shot = 88741;
% shot = 88742; 
% shot = 88647

shot = 88903

% possible shots (4/12/2025)
shots = [88564  88565  88566  88567  88577  88579  88582  88583  88584  88586  88646  88647  88734  88736  88737  88738  88739  88740  88741  88742];

mdsconnect('tcvdata.epfl.ch');
successful_shots = 0;
for shot = shots
    fprintf('Running analysis for shot %d...\n', shot);
    try

    % START_TIME = 0.05 % start time for analysis
    START_TIME = 0.0 %0.2 % start time for analysis
    STOP_TIME = 3.0 %1.6 %3.0; % stop time for analysis
    dec = 10; % decimation factor 10

    
    nrt = get_scd_mems(shot,'liuqenn');
    Fxn = nrt.liuqenn_Fxn;
    Brn = nrt.liuqenn_Brn;
    Bzn = nrt.liuqenn_Bzn;
    Fx = Fxn.Data; Br = Brn.Data; Bz = Bzn.Data;
    
    % figure;
    % subplot(3,1,1); grid on; plot(Fxn); title('Fxn');
    % subplot(3,1,2); grid on; plot(Brn); title('Brn');
    % subplot(3,1,3); grid on; plot(Bzn); title('Bzn');
    
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
    % swap 17,18,19 to correct order
    Ia_17 = Ia(:,17);
    Ia_18 = Ia(:,18);
    Ia_19 = Ia(:,19); 
    Ia(:,17) = Ia_19;
    Ia(:,18) = Ia_17;
    Ia(:,19) = Ia_18;
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
    
    Fx_or = Fx; Br_or = Br; Bz_or = Bz; % save original online outputs before decimation for later comparison
    t_or = t; phys_or = phys;
    
    
    %% analysis
    
    
    % decimate
    Fx = Fx(1:dec:end,:);
    Br = Br(1:dec:end,:);
    Bz = Bz(1:dec:end,:);
    t = t(1:dec:end);
    phys = phys(1:dec:end, :);
    
    phys_RT = phys;
    
    rt_lx.Bm = phys_RT(:,1:38)';
    rt_lx.Ff = phys_RT(:,39:76)';
    rt_lx.Ia = phys_RT(:,77:95)';
    rt_lx.Uf = phys_RT(:,96:133)';
    rt_lx.rBt = phys_RT(:,134)';
    
    % INTERP_METHOD = 'linear'; % interpolation method for interp2
    % INTERP_METHOD = 'spline'; % interpolation method for interp2
    INTERP_METHOD = 'cubic'; % interpolation method for interp2
    
    r = [1.0864 1.0313 0.9357 0.8253 0.7297 0.6746 0.6746 0.7297 0.8253 0.9357 1.0313 1.0864 1.0137 0.9747 0.8805 0.7863 0.7473 0.7863 0.8805 0.9747 0.9182 0.9182 0.8428 0.8428 0.8805];
    z = [0.1574  0.4299  0.5873  0.5873  0.4299  0.1574 -0.1574 -0.4299 -0.5873 -0.5873 -0.4299 -0.1574 0. 0.2687 0.38 0.2687 0. -0.2687 -0.38 -0.2687 -0.1075 0.1075 0.1075 -0.1075  0.];
    
    % online inference captured in real time
    Fx_online = Fx; Br_online = Br; Bz_online = Bz;
    
    fprintf('Running LIUQE versions...\n');
    
    [L_map, LY_map] = mds2meq(shot, 'LIUQE.M'); % full maps from the precalculated LIUQE output
    Fx_map = LY_map.Fx; Iy_map = LY_map.Iy; t_map = LY_map.t; Ip_map = LY_map.Ip;
    r_map = L_map.rx; z_map = L_map.zx;
    assert(size(Fx_map, 3) == numel(t_map) && size(Iy_map, 3) == numel(t_map) && numel(Ip_map) == numel(t_map), 'Precalculated LIUQE output sizes do not match');

    [L,LX] = liuqe(shot,t);
    LX.Bm=rt_lx.Bm; LX.Ff=rt_lx.Ff; LX.Uf=rt_lx.Uf; LX.Ia=rt_lx.Ia;% LX.rBt=rt_lx.rBt; % use phys_RT inputs
    LX_liuqe = LX; % save for later
    
    % rtliuqe (approx)
    % [Lrt,LXrt,LYrt] = liuqe(shot,t,'itert',1);
    % LYrt = liut(L,LX, 'itert', 1); % double check that itert=1 is the same as rtliuqe (too good like this)
    % [Lrt,LXrt,LYrt] = liuqe(shot,t, 'LX', LX, 'itert',1); % not 100% correct
    % Lrt = liuqe(shot, [], 'itert', 1); LYrt = liut(Lrt, LX, 'itert', 1); % better version

    % MATLAB emulation of RT-LIUQE using the parameters deployed for this shot.
    [Lrt,LXrt] = liuqe(shot,t, 'psrc','RTLIUQE', 'slx',1, 'bslvdprec',true);
    assert(Lrt.P.itert > 0, 'RT-LIUQE parameters did not enable itert');
    assert(Lrt.P.ivesm && startsWith(Lrt.P.selu, 's'), 'Expected RT-LIUQE to use measured segment currents');
    LXrt.Bm = rt_lx.Bm;
    LXrt.Ff = rt_lx.Ff;
    LXrt.Ia = rt_lx.Ia;
    LXrt.Uf = rt_lx.Uf;
    LXrt.Iu = LXrt.Uf ./ (-Lrt.G.Ru(:));
    Is = LXrt.Uf ./ (-Lrt.G.Rs(:));
    LXrt.Ip = Lrt.G.Ipm*LXrt.Bm + Lrt.Ipa*LXrt.Ia + Lrt.Ips*Is;
    % Keep the measured LXrt.rBt instead of replacing it with the artificial zeros in rt_lx.rBt.
    LXrt = liux(Lrt,LXrt);

    Lrt.liurtemu = true;
    LYrt = liut(Lrt,LXrt, 'rst',true, 'slx',false);
    [Fx_lrt, Br_lrt, Bz_lrt] = interp_liuqe(r,z, Lrt, LXrt, LYrt, INTERP_METHOD);
    
    % lih
    [Lih,LXih,LYih] = lih('tcv', shot, t, 'debug', 1);
    [Fx_lih, Br_lih, Bz_lih] = interp_liuqe(r,z, L, LX, LYih, INTERP_METHOD);
    
    % standard liuqe
    LY = liut(L,LX);
    [Fx_liu, Br_liu, Bz_liu] = interp_liuqe(r,z, L, LX, LY, INTERP_METHOD);
    
    % small and turbo recalculated on phys_RT inputs
    [Fx_small2, Br_small2, Bz_small2] = liuqe_nn_LX(rt_lx, single(r'), single(z'));
    [Fx_turbo2, Br_turbo2, Bz_turbo2] = liuqe_nn_LX(rt_lx, single(r'), single(z'), 'turbo');
    
    % small and turbo recalculated on LX
    [Fx_small, Br_small, Bz_small] = liuqe_nn_LX(LX_liuqe, single(r'), single(z')); % this assigns to base phys_LX
    [Fx_turbo, Br_turbo, Bz_turbo] = liuqe_nn_LX(LX_liuqe, single(r'), single(z'), 'turbo');
    
    
    % print abs err, std of the err and max err of Fx, Br, Bz between liuqe and liuqenn
    e_Fx_online = abs(Fx_liu - Fx_online);
    e_Br_online = abs(Br_liu - Br_online);
    e_Bz_online = abs(Bz_liu - Bz_online);
    fprintf('LIUQENN ONLINE vs LIUQE:\n');
    fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_online(:)), std(e_Fx_online(:)), max(e_Fx_online(:)));
    fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_online(:)), std(e_Br_online(:)), max(e_Br_online(:)));
    fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_online(:)), std(e_Bz_online(:)), max(e_Bz_online(:)));
    
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
    
    % small on LX vs liuqe
    e_Fx_small = abs(Fx_small - Fx_liu);
    e_Br_small = abs(Br_small - Br_liu);
    e_Bz_small = abs(Bz_small - Bz_liu);
    fprintf('SMALL vs LIUQE:\n');
    fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_small(:)), std(e_Fx_small(:)), max(e_Fx_small(:)));
    fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_small(:)), std(e_Br_small(:)), max(e_Br_small(:)));
    fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_small(:)), std(e_Bz_small(:)), max(e_Bz_small(:)));

    % small on phys_RT vs liuqe
    e_Fx_small2 = abs(Fx_small2 - Fx_liu);
    e_Br_small2 = abs(Br_small2 - Br_liu);
    e_Bz_small2 = abs(Bz_small2 - Bz_liu);
    fprintf('SMALL on phys_RT vs LIUQE:\n');
    fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_small2(:)), std(e_Fx_small2(:)), max(e_Fx_small2(:)));
    fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_small2(:)), std(e_Br_small2(:)), max(e_Br_small2(:)));
    fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_small2(:)), std(e_Bz_small2(:)), max(e_Bz_small2(:)));

    % turbo on LX vs liuqe
    e_Fx_turbo = abs(Fx_turbo - Fx_liu);
    e_Br_turbo = abs(Br_turbo - Br_liu);
    e_Bz_turbo = abs(Bz_turbo - Bz_liu);
    fprintf('TURBO vs LIUQE:\n');
    fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_turbo(:)), std(e_Fx_turbo(:)), max(e_Fx_turbo(:)));
    fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_turbo(:)), std(e_Br_turbo(:)), max(e_Br_turbo(:)));
    fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_turbo(:)), std(e_Bz_turbo(:)), max(e_Bz_turbo(:)));

    % turbo on phys_RT vs liuqe
    e_Fx_turbo2 = abs(Fx_turbo2 - Fx_liu);
    e_Br_turbo2 = abs(Br_turbo2 - Br_liu);
    e_Bz_turbo2 = abs(Bz_turbo2 - Bz_liu);
    fprintf('TURBO on phys_RT vs LIUQE:\n');
    fprintf('Fx: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Fx_turbo2(:)), std(e_Fx_turbo2(:)), max(e_Fx_turbo2(:)));
    fprintf('Br: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Br_turbo2(:)), std(e_Br_turbo2(:)), max(e_Br_turbo2(:)));
    fprintf('Bz: mean err: %.2e, std err: %.2e, max err: %.2e\n', mean(e_Bz_turbo2(:)), std(e_Bz_turbo2(:)), max(e_Bz_turbo2(:)));

    % save every numeric array as single precision
    save_names = {'t_or', 'phys_or', 'Fx_or', 'Br_or', 'Bz_or', 'Fx_online', 'Br_online', 'Bz_online', ...
        'Fx_lrt', 'Br_lrt', 'Bz_lrt', 'Fx_lih', 'Br_lih', 'Bz_lih', 'Fx_liu', 'Br_liu', 'Bz_liu', ...
        'Fx_small', 'Br_small', 'Bz_small', 'Fx_small2', 'Br_small2', 'Bz_small2', ...
        'Fx_turbo', 'Br_turbo', 'Bz_turbo', 'Fx_turbo2', 'Br_turbo2', 'Bz_turbo2', 't', 'phys', ...
        'Fx_map', 'Iy_map', 't_map', 'Ip_map', 'r_map', 'z_map'};
    save_data = struct();
    for i = 1:numel(save_names)
        save_data.(save_names{i}) = single(eval(save_names{i}));
    end
    save(sprintf('realtime_analysis/rt_shot%d.mat', shot), '-struct', 'save_data');
    successful_shots = successful_shots + 1;
            
    catch ME
        fprintf(2, 'Skipping shot %d: %s\n', shot, ME.message);
        continue;
    end
end
mdsdisconnect;

if successful_shots == 0
    fprintf('No shots completed successfully; skipping plots.\n');
    return;
end

%plotting
plt_dec = 2; % number of points to plot

figure;
subplot(3,1,1); grid on; hold on;
plot(t, Fx_online(:,1:plt_dec:end));
% plot(t, Fx_lrt(:,1:plt_dec:end), '--');
plot(t, Fx_liu(:,1:plt_dec:end), ':');
title('Fx'); xlabel('Time (s)');
subplot(3,1,2); grid on; hold on;
plot(t, Br_online(:,1:plt_dec:end));
% plot(t, Br_lrt(:,1:plt_dec:end), '--');
plot(t, Br_liu(:,1:plt_dec:end), ':');
title('Br'); xlabel('Time (s)');
subplot(3,1,3); grid on; hold on;
plot(t, Bz_online(:,1:plt_dec:end));
% plot(t, Bz_lrt(:,1:plt_dec:end), '--');
plot(t, Bz_liu(:,1:plt_dec:end), ':');
title('Bz'); xlabel('Time (s)');

% comparison of the phys inputs
figure;
bm_rt = phys_RT(:,1:38); ff_rt = phys_RT(:,39:76); ia_rt = phys_RT(:,77:95); uf_rt = phys_RT(:,96:133); rbt_rt = phys_RT(:,134);
bm_lx = phys_LX(:,1:38); ff_lx = phys_LX(:,39:76); ia_lx = phys_LX(:,77:95); uf_lx = phys_LX(:,96:133); rbt_lx = phys_LX(:,134);
e_bm=bm_rt-bm_lx; e_ff=ff_rt-ff_lx; e_ia=ia_rt-ia_lx; e_uf=uf_rt-uf_lx; % pure error
% e_bm=abs(e_bm); e_ff=abs(e_ff); e_ia=abs(e_ia); e_uf=abs(e_uf); % absolute error

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
    % Same as for Brs'; %
    Bz(:,end    ,:) =  i4pirdr(end)     .* (+Fx(:,end-2,:) - 4*Fx(:,end-1,:) + 3*Fx(:,end,:));
    Bz(:,1      ,:) =  i4pirdr(1)       .* (-Fx(:,    3,:) + 4*Fx(:,    2,:) - 3*Fx(:,  1,:));
end
