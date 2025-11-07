function [d] = ds_prep_save_shot(shot, LX, LY, t_mask, dec, save_dir)

    %% extract quantities
    % Ouputs
    Fx = LY.Fx; % Plasma poloidal flux map | `(rx,zx,t)` | `[Wb]` |
    Iy = LY.Iy; % Plasma current density map | `(ry,zy,t)` | `[A/m^2]` |

    % Inputs 
    Bm = LX.Bm; 
    Ff = LX.Ff;
    Ia = LX.Ia; 
    Ip = LX.Ip;
    Uf = LX.Uf;
    rBt = LX.rBt; 

    % check the dimensions
    assert(all(size(LX.t) == size(LY.t)), 'Time vector sizes do not match');
    t = LX.t;
    nt = numel(t); % number of time samples
    assert(all(size(Fx) == [65, 28, nt]), 'Fx has wrong size');
    assert(all(size(Iy) == [63, 26, nt]), 'Iy has wrong size');
    % filter out the NaN/Inf values [MILD]
    mFx = reshape(all(all(~isnan(Fx) & ~isinf(Fx),1),2), [], 1);
    mIy = reshape(all(all(~isnan(Iy) & ~isinf(Iy),1),2), [], 1);

    mBm = reshape(all(~isnan(Bm) & ~isinf(Bm), 1), [], 1);
    mFf = reshape(all(~isnan(Ff) & ~isinf(Ff), 1), [], 1);
    mIa = reshape(all(~isnan(Ia) & ~isinf(Ia), 1), [], 1);
    mIp = reshape(~isnan(Ip) & ~isinf(Ip), [], 1);
    mUf = reshape(all(~isnan(Uf) & ~isinf(Uf), 1), [], 1);
    mrBt = reshape(~isnan(rBt) & ~isinf(rBt), [], 1);
    t_mask = reshape(t_mask, [], 1);

    valid = mFx & mIy & mBm & mFf & mIa & mIp & mUf & mrBt & t_mask; % keep only the samples where all quantities are valid + additional time mask
    fprintf('\tvalid samples -> [TOT: %.1f%%, %d] Fx:%.1f%%, Iy: %.1f%%, \n\tBm: %.1f%%, Ff: %.1f%%, Ia: %.1f%%, Ip: %.1f%%, Uf: %.1f%%, rBt: %.1f%%\n', ...
        100*sum(valid)/nt, sum(valid), 100*sum(mFx)/nt, 100*sum(mIy)/nt, 100*sum(mBm)/nt, 100*sum(mFf)/nt, 100*sum(mIa)/nt, 100*sum(mIp)/nt, 100*sum(mUf)/nt, 100*sum(mrBt)/nt);
    assert(sum(valid) > 0.5 * nt, 'Nan/Inf filter -> not enough valid samples');
    
    % keep only the valid samples
    t = t(valid);
    Fx = Fx(:,:,valid);
    Iy = Iy(:,:,valid);
    Bm = Bm(:,valid);
    Ff = Ff(:,valid);
    Ia = Ia(:,valid);
    Ip = Ip(valid);
    Uf = Uf(:,valid);
    rBt = rBt(valid);

    % decimate
    t = t(1:dec:end);
    Fx = Fx(:,:,1:dec:end);
    Iy = Iy(:,:,1:dec:end);
    Bm = Bm(:,1:dec:end);
    Ff = Ff(:,1:dec:end);
    Ia = Ia(:,1:dec:end);
    Ip = Ip(1:dec:end);
    Uf = Uf(:,1:dec:end);
    rBt = rBt(1:dec:end);

    % convert to single precision to save space
    Fx = single(Fx); 
    Iy = single(Iy);
    Bm = single(Bm);
    Ff = single(Ff);
    Ia = single(Ia);
    Ip = single(Ip);
    Uf = single(Uf);
    rBt = single(rBt);

    % print final sizes
    fprintf('\tFinal sizes -> Fx:%s, Iy:%s, \n\tBm:%s, Ff:%s, Ia:%s, Ip:%s, Uf:%s, rBt:%s\n', ...
        mat2str(size(Fx)), mat2str(size(Iy)), mat2str(size(Bm)), mat2str(size(Ff)), mat2str(size(Ia)), mat2str(size(Ip)), mat2str(size(Uf)), mat2str(size(rBt)));

    % save data into a .mat file
    save_file = fullfile(save_dir, sprintf('%d.mat', shot));

    save(save_file, 't', 'Fx', 'Iy', 'Bm', 'Ff', 'Ia', 'Ip', 'Uf', 'rBt');
    
    d.t=t; d.Fx=Fx; d.Iy=Iy; d.Bm=Bm; d.Ff=Ff; d.Ia=Ia; d.Ip=Ip; d.Uf=Uf; d.rBt=rBt;
    
    fprintf('\x1b[32m\tData saved to: %s\x1b[0m\n', save_file);
end