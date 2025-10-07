clear all; close all; clc;

% TODO consider adding vessel currents: I_VESSEL

START_SHOT = 77662; % Dec 2022, https://spcwiki.epfl.ch/wiki/Alma_database
% END_SHOT = 85804; % April 2025
END_SHOT = 87899; % Sept 2025
% N_SHOTS = 100; % Number of shots to process
N_SHOTS = END_SHOT-START_SHOT; % Number of shots to process

% Directory to save the output .mat files
% OUT_DIR = 'ds'; % testing
OUT_DIR = '/NoTivoli/grandin/ds'; % more space available

DECIMATION = 10; % 4 Decimation factor for the time vector

MIN_TIME_SAMPLES = 5; % 10 Minimum number of time samples to keep the shot
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

shots = randi([START_SHOT, END_SHOT], 1, N_SHOTS);
total_shots = 0;
fprintf('Shots: %s\n', mat2str(shots));
fprintf('\nStarting data retrieval loop...\n');
shot_processing_times = zeros(1, length(shots));
for i = 1:length(shots)
    start_time = tic; % Start timer for processing each shot

    shot = shots(i);
    fprintf('\x1b[33mProcessing shot %d [%d/%d]\x1b[0m\n', shot, i, length(shots));

    try 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load liuqe data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t = liuqe_times(shot); % get times
        t = t(1:DECIMATION:end); % decimate times
        
        % run liuqe
        [L, LX] = liuqe(shot,t,'LYall',false,'ndom',4);
        
        assert(all(abs(LX.t - t) < 1e-6), 'LX.t is not sufficiently close to t');
        
        LY = liut(L, LX);
        
        valid_t = ismember(LX.t, LY.t); % match times
        t = LX.t(valid_t); % keep only valid times
        nt = length(t);
        assert(nt >= MIN_TIME_SAMPLES, 'Not enough time samples');
        fprintf('\tNumber of times: %d\n', nt);

        %% extract quantities
        % Ouputs
        Fx = LY.Fx; % Plasma poloidal flux map | `(rx,zx,t)` | `[Wb]` |
        Iy = LY.Iy; % Plasma current density map | `(ry,zy,t)` | `[A/m^2]` |

        % Inputs
        Bm = LX.Bm(:,valid_t); % magnetic probe measurements | (*,t) | [T]
        Ff = LX.Ff(:,valid_t); % flux loop poloidal flux | (*,t) | [Wb] 
        Ia = LX.Ia(:,valid_t); % poloidal field coil currents | (*,t) | [A]
        Ip = LX.Ip(valid_t);   % plasma current | (t) | [A]

        % check the dimensions and show actual size if wrong
        assert(all(size(Fx) == [65, 28, nt]), sprintf('Fx has wrong size: %s', mat2str(size(Fx))));
        assert(all(size(Iy) == [63, 26, nt]), sprintf('Iy has wrong size: %s', mat2str(size(Iy))));

        assert(all(size(Bm) == [38, nt]), sprintf('Bm has wrong size: %s', mat2str(size(Bm))));
        assert(all(size(Ff) == [38, nt]), sprintf('Ff has wrong size: %s', mat2str(size(Ff))));
        assert(all(size(Ia) == [19, nt]), sprintf('Ia has wrong size: %s', mat2str(size(Ia))));
        assert(all(size(Ip) == [1, nt]), sprintf('Ip has wrong size: %s', mat2str(size(Ip))));

        % filter out the NaN/Inf values [MILD]
        mFx = reshape(all(all(~isnan(Fx) & ~isinf(Fx),1),2), [], 1);
        mIy = reshape(all(all(~isnan(Iy) & ~isinf(Iy),1),2), [], 1);

        mBm = reshape(all(~isnan(Bm) & ~isinf(Bm), 1), [], 1);
        mFf = reshape(all(~isnan(Ff) & ~isinf(Ff), 1), [], 1);
        mIa = reshape(all(~isnan(Ia) & ~isinf(Ia), 1), [], 1);
        mIp = reshape(~isnan(Ip) & ~isinf(Ip), [], 1);

        % valid = mFx & mIy & mBr & mBz & ...
        valid = mFx & mIy & mBm & mFf & mIa & mIp;
        fprintf('\tvalid samples -> [TOT: %.1f%%, %d] Fx:%.1f%%, Iy: %.1f%%, \n\tBm: %.1f%%, Ff: %.1f%%, Ia: %.1f%%, Ip: %.1f%%\n', ...
            100*sum(valid)/nt, sum(valid), 100*sum(mFx)/nt, 100*sum(mIy)/nt, ...
            100*sum(mBm)/nt, 100*sum(mFf)/nt, 100*sum(mIa)/nt, 100*sum(mIp)/nt);
        assert(sum(valid) > 0.5 * nt, 'Nan/Inf filter -> not enough valid samples');
        
        % keep only the valid samples
        t = t(valid);
        Fx = Fx(:,:,valid);
        Iy = Iy(:,:,valid);
        Bm = Bm(:,valid);
        Ff = Ff(:,valid);
        Ia = Ia(:,valid);
        Ip = Ip(valid);

        % convert to single precision to save space
        Fx = single(Fx); 
        Iy = single(Iy);
        Bm = single(Bm);
        Ff = single(Ff);
        Ia = single(Ia);
        Ip = single(Ip);

        % print final sizes
        fprintf('\tFinal sizes -> Fx:%s, Iy:%s, \n\tBm:%s, Ff:%s, Ia:%s, Ip:%s\n', ...
            mat2str(size(Fx)), mat2str(size(Iy)), ...
            mat2str(size(Bm)), mat2str(size(Ff)), mat2str(size(Ia)), mat2str(size(Ip)));
        mdsclose; % Close the MDSplus connection
        total_shots = total_shots + 1;

        % save data into a .mat file
        save_file = fullfile(OUT_DIR, sprintf('%d.mat', shot));
        save(save_file, 't', 'Fx', 'Iy', 'Bm', 'Ff', 'Ia', 'Ip');
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