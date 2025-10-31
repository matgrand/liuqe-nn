clear all; close all; clc;

START_SHOT = 77662; % Dec 2022, https://spcwiki.epfl.ch/wiki/Alma_database
% END_SHOT = 85804; % April 2025
END_SHOT = 87899; % Sept 2025
% N_SHOTS = 100; % Number of shots to process
N_SHOTS = END_SHOT-START_SHOT; % Number of shots to process
N_SHOTS = 5; % Number of shots to process

% Directory to save the output .mat files
OUT_DIR = '/NoTivoli/grandin/ds'; % more space available on NoTivoli

DECIMATION = 6; % 6 Decimation factor for the time vector

MIN_TIME_SAMPLES = 10; % Minimum number of time samples to keep the shot
MAX_IP_PERC_DIFF = 2.5; % Maximum percentage difference between IPLIUQE and IP

fprintf('Parameters:\n');
fprintf('\tSTART_SHOT: %d\n', START_SHOT);
fprintf('\tEND_SHOT: %d\n', END_SHOT);
fprintf('\tN_SHOTS: %d\n', N_SHOTS);
fprintf('\tOUT_DIR: %s\n', OUT_DIR);
fprintf('\tDECIMATION: %d\n', DECIMATION);
fprintf('\tMIN_TIME_SAMPLES: %d\n', MIN_TIME_SAMPLES);
fprintf('\tMAX_IP_PERC_DIFF: %.2f\n', MAX_IP_PERC_DIFF);
fprintf('\n');

fprintf('Deleting old figures...\n');
delete(fullfile('figs', 'ip_*.png')); % Delete old figures

if ~exist(OUT_DIR, 'dir') mkdir(OUT_DIR); fprintf('Output directory created: %s\n', OUT_DIR);
else delete(fullfile(OUT_DIR, '*')); fprintf('Output directory already exists. Old files deleted: %s\n', OUT_DIR);
end % Create output directory if it doesn't exist

mdsconnect('tcvdata.epfl.ch'); % Connect to the MDSplus server

shots = randi([START_SHOT, END_SHOT], 1, N_SHOTS);
total_shots = 0;
fprintf('Shots: %s\n', mat2str(shots));
fprintf('\nStarting data retrieval loop...\n');
for i = 1:length(shots)
    start_time = tic; % Start timer for processing each shot

    shot = shots(i);
    fprintf('\x1b[33mProcessing shot %d (%d of %d)\x1b[0m\n', shot, i, length(shots));

    try 
        mdsopen('tcv_shot', shot); % Open the MDSplus connection to the TCV database

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load liuqe data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [L, LY] = mds2meq(shot, 'LIUQE.M'); % get liuqe outputs from mdsplus
        [L, LX] = liuqe(shot, LY.t); % get liuqe inputs 
        
        t = LY.t'; % time vector
        [t2, ip2] = tcvget('IP', t); % calculated using magnetics at liuqe times

        % analyze the time vector
        assert(max(abs(t2 - t)) < 1e-8, 'Time vectors do not coincide');
	    assert(numel(t) > 1, sprintf('Time vector has insufficient elements: t:%s', mat2str(size(t))));
        nt = numel(t); % number of time samples
        t_diff = abs(t - t2);
        fprintf('\ttime steps -> n: %d, mean: %.2f [µs], std: %.2f [µs]\n', numel(t), mean(diff(t) * 1e6), std(diff(t) * 1e6));
        assert(max(t_diff) < 1e-8, 'Times do not coincide');
        
        % analyze the plasma current
        ip1 = LY.Ip'; % IPLIUQE
	    assert(all(size(ip1) == size(ip2)), 'Ip sizes wrong');
        avg_ip = mean(abs(ip2));
        ip_diff = abs(ip1 - ip2);
        avg_diff = mean(ip_diff);
        perc_diff = ip_diff ./ abs(ip2) * 100;
        fprintf('\tip average difference -> %.2f [A] (%.1f%%)\n', avg_diff, mean(perc_diff));
        % keep only the samples where IP and IPLIUQE are similar
        % assert(mean(perc_diff) < MAX_IP_PERC_DIFF, 'Difference between IPLIUQE and IP is too high'); % very strict
        ip_valid1 = perc_diff < MAX_IP_PERC_DIFF;
        fprintf('\tip filtered -> %.1f%%, remaining -> %d/%d \n', 100*(1-sum(ip_valid1)/numel(ip_valid1)), sum(ip_valid1), numel(ip_valid1));
        assert(sum(ip_valid1) > 0.8 * numel(ip_valid1), 'IP MEAS and IPLIUQE are different in too many samples');

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

        valid = mFx & mIy & mBm & mFf & mIa & mIp & mUf & mrBt & ip_valid1; % keep only the samples where all quantities are valid and IP is valid
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
        t = t(1:DECIMATION:end);
        Fx = Fx(:,:,1:DECIMATION:end);
        Iy = Iy(:,:,1:DECIMATION:end);
        Bm = Bm(:,1:DECIMATION:end);
        Ff = Ff(:,1:DECIMATION:end);
        Ia = Ia(:,1:DECIMATION:end);
        Ip = Ip(1:DECIMATION:end);
        Uf = Uf(:,1:DECIMATION:end);
        rBt = rBt(1:DECIMATION:end);

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
        mdsclose; % Close the MDSplus connection
        total_shots = total_shots + 1;

        % save data into a .mat file
        save_file = fullfile(OUT_DIR, sprintf('%d.mat', shot));

        save(save_file, 't', 'Fx', 'Iy', 'Bm', 'Ff', 'Ia', 'Ip', 'Uf', 'rBt');
        fprintf('\x1b[32m\tData saved to: %s\x1b[0m\n', save_file);
        
        shot_processing_times(i) = toc(start_time);
        fprintf('\tProc time: %.2f s, ETA: %.0f min\n', shot_processing_times(i), sum(shot_processing_times) / i * (length(shots) - i) / 60);
        
    catch ME
        fprintf('\x1b[31m\tError processing shot %d: %s\x1b[0m\n', shot, ME.message);
        % for k = 1:length(ME.stack)
        %     fprintf('\t\tIn %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
        % end
        continue; % Skip to the next shot on error
    end % try-catch
end % end shots loop

mdsdisconnect; % Disconnect from MDSplus

fprintf('\nProcessing complete for all shots.\n');
fprintf('Output files saved in: %s\n', OUT_DIR);

% save a random mat in the /NoTivoli/grandin/ directory to show it is finished
save(fullfile('/NoTivoli/grandin/', 'ds_done.mat'), 'shot_processing_times');
