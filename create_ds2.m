clear all; close all; clc;

% TODO consider adding vessel currents: I_VESSEL

START_SHOT = 77662; % Dec 2022, https://spcwiki.epfl.ch/wiki/Alma_database
% END_SHOT = 85804; % April 2025
% START_SHOT = 87700; % May 2025
END_SHOT = 87899; % Sept 2025
% N_SHOTS = 10; % Number of shots to process
N_SHOTS = END_SHOT-START_SHOT; % Number of shots to process

% Directory to save the output .mat files
% OUT_DIR = 'ds'; % testing
OUT_DIR = '/NoTivoli/grandin/ds'; % more space available

DECIMATION = 40; % 4 Decimation factor for the time vector
NDOM = 4; % Number of domains for liuqe (must be >=3 for doublet)
% ITERA = 50; % Max Number of iterations for liuqe NOTE: cant change this bc has some weird values for specific shapes (like 888 for doublets)
ITERFRZ = 20; %20; % Number of Picard iterations before which plasma domain is frozen to ease convergence (inactive if itert>0)
ITERH = 20; %20; % Max Number of FE solver iterations
ITERQ = 10; % iterations for the LCFS

MIN_TIME_SAMPLES = 5; % 10 Minimum number of time samples to keep the shot

fprintf('Parameters:\n');
fprintf('\tSTART_SHOT: %d\n', START_SHOT);
fprintf('\tEND_SHOT: %d\n', END_SHOT);
fprintf('\tN_SHOTS: %d\n', N_SHOTS);
fprintf('\tOUT_DIR: %s\n', OUT_DIR);
fprintf('\tDECIMATION: %d\n', DECIMATION);
fprintf('\tNDOM: %d\n', NDOM);
% fprintf('\tITERA: %d\n', ITERA);
fprintf('\tITERFRZ: %d\n', ITERFRZ);
fprintf('\tITERH: %d\n', ITERH);
fprintf('\tITERQ: %d\n', ITERQ);
fprintf('\tMIN_TIME_SAMPLES: %d\n', MIN_TIME_SAMPLES);
fprintf('\n');

fprintf('Deleting old figures...\n');
delete(fullfile('figs', 'ip_*.png')); % Delete old figures

if ~exist(OUT_DIR, 'dir') mkdir(OUT_DIR); fprintf('Output directory created: %s\n', OUT_DIR);
else delete(fullfile(OUT_DIR, '*')); fprintf('Output directory already exists. Old files deleted: %s\n', OUT_DIR);
end % Create output directory if it doesn't exist

shots = randi([START_SHOT, END_SHOT], 1, N_SHOTS);
% shots = [87835, 87851];
% shots = [87835];
% shots = [87851];
total_shots = 0;
fprintf('Shots: %s\n', mat2str(shots));
fprintf('\nStarting data retrieval loop...\n');
shot_processing_times = zeros(1, length(shots));
for i = 1:length(shots)
    start_time = tic; % Start timer for processing each shot

    shot = shots(i);
    fprintf('\x1b[33mProcessing shot %d [%d/%d]\x1b[0m\n', shot, i, length(shots));

    try 
        %% load liuqe data
        t = liuqe_times(shot); % get times
        t = t(1:DECIMATION:end); % decimate times
        
        %% run liuqe
        [L, LX, LY] = liuqe(shot,t,'ndom',NDOM,'iterh',ITERH,'iterfrz',ITERFRZ,'iterq',ITERQ);

        t = LX.t;
        nt = length(t);
        assert(nt >= MIN_TIME_SAMPLES, 'Not enough time samples');

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
        
        % check the dimensions
        assert(all(size(Fx) == [65, 28, nt]), sprintf('Fx has wrong size: %s', mat2str(size(Fx))));
        assert(all(size(Iy) == [63, 26, nt]), sprintf('Iy has wrong size: %s', mat2str(size(Iy))));
        assert(all(size(Bm) == [38, nt]), sprintf('Bm has wrong size: %s', mat2str(size(Bm))));
        assert(all(size(Ff) == [38, nt]), sprintf('Ff has wrong size: %s', mat2str(size(Ff))));
        assert(all(size(Ia) == [19, nt]), sprintf('Ia has wrong size: %s', mat2str(size(Ia))));
        assert(all(size(Ip) == [1, nt]), sprintf('Ip has wrong size: %s', mat2str(size(Ip))));
        assert(all(size(rq) == [32*NDOM, nt]), sprintf('rq has wrong size: %s', mat2str(size(rq))));
        assert(all(size(zq) == [32*NDOM, nt]), sprintf('zq has wrong size: %s', mat2str(size(zq))));

        % filter out the NaN/Inf values
        assert(~any(isnan(Fx(:))) && ~any(isinf(Fx(:))), 'Fx contains NaN or Inf values');
        assert(~any(isnan(Iy(:))) && ~any(isinf(Iy(:))), 'Iy contains NaN or Inf values');
        assert(~any(isnan(Bm(:))) && ~any(isinf(Bm(:))), 'Bm contains NaN or Inf values');
        assert(~any(isnan(Ff(:))) && ~any(isinf(Ff(:))), 'Ff contains NaN or Inf values');
        assert(~any(isnan(Ia(:))) && ~any(isinf(Ia(:))), 'Ia contains NaN or Inf values');
        assert(~any(isnan(Ip(:))) && ~any(isinf(Ip(:))), 'Ip contains NaN or Inf values');
        assert(~any(isnan(rq(:))) && ~any(isinf(rq(:))), 'rq contains NaN or Inf values');
        assert(~any(isnan(zq(:))) && ~any(isinf(zq(:))), 'zq contains NaN or Inf values');

        % print final sizes
        fprintf('\tFinal sizes -> Fx:%s, Iy:%s, Bm:%s, Ff:%s, Ia:%s, Ip:%s, rq:%s, zq:%s\n', ...
            mat2str(size(Fx)), mat2str(size(Iy)), mat2str(size(Bm)), mat2str(size(Ff)), mat2str(size(Ia)), mat2str(size(Ip)), mat2str(size(rq)), mat2str(size(zq)));
        total_shots = total_shots + 1;

        % save data into a .mat file
        save_file = fullfile(OUT_DIR, sprintf('%d.mat', shot));
        save(save_file, 't', 'Fx', 'Iy', 'Bm', 'Ff', 'Ia', 'Ip', 'rq', 'zq');
        fprintf('\x1b[32m\tData saved to: %s\x1b[0m\n', save_file);
        
        shot_processing_times(i) = toc(start_time);
        fprintf('\tProc time: %.2f s, ETA: %.0f min\n', shot_processing_times(i), sum(shot_processing_times) / i * (length(shots) - i) / 60);
        
        % % Quick plot of Fx and Iy at mid time
        % figure; mid = round(nt/2); sgtitle(sprintf('Shot %d', shot));
        % [rr, zz] = meshgrid(L.rx, L.zx);
        % subplot(1,2,1); contourf(rr, zz, Fx(:,:,mid), 20, 'LineColor', 'none'); hold on;
        % scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled'); axis xy; colorbar; title('Fx'); 
        % [rry, zzy] = meshgrid(L.ry, L.zy);
        % subplot(1,2,2); contourf(rry, zzy, Iy(:,:,mid), 20, 'LineColor', 'none'); hold on;
        % scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled'); axis xy; colorbar; title('Iy');


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