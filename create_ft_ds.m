clear all; close all; clc;

% set warning off for everything
warning('off', 'all');

% script to create a fine-tuning dataset 

shots = [87702 87703 87704 87706 87707 87708 87709 87710 87711 87712 87713 87714 87714 87715 87716 87717 87718 87719 87732 87735 87736 87737 87738 87739 87740 87741 87742 87743 87744 87745 87746 87747 87748 87750 87751 87752 87753 87755 87767 87767 87768 87769 87770 87771 87772 87773 87774 87775 87776 87777 87778 87779 87793 87796 87796 87797 87798 87799 87800 87801 87802 87803 87804 87805 87806 87807 87808 87809 87810 87811 87812 87828 87829 87832 87833 87834 87835 87836 87837 87838 87839 87840 87841 87842 87843 87844 87845 87846 87848 87849 87850 87866 87867 87868 87869 87870 87871 87872 87873 87874 87875 87876 87877 87878 87879 87880 87881 87882 87882 87892 87893 87895 87896 87897 87898 87899];

% Directory to save the output .mat files
OUT_DIR = '/NoTivoli/grandin/ft_ds'; % more space available

DECIMATION = 3; %2; % 4 Decimation factor for the time vector
NDOM = 3; %4; % Number of domains for liuqe (must be >=3 for doublet)
ITERFRZ = 20; %20; % Number of Picard iterations before which plasma domain is frozen to ease convergence (inactive if itert>0)
ITERH = 20; %20; % Max Number of FE solver iterations

MIN_TIME_SAMPLES = 5; % 10 Minimum number of time samples to keep the shot

fprintf('Parameters:\n');
fprintf('\tOUT_DIR: %s\n', OUT_DIR);
fprintf('\tDECIMATION: %d\n', DECIMATION);
fprintf('\tNDOM: %d\n', NDOM);
fprintf('\tITERFRZ: %d\n', ITERFRZ);
fprintf('\tITERH: %d\n', ITERH);
fprintf('\tMIN_TIME_SAMPLES: %d\n', MIN_TIME_SAMPLES);
fprintf('\n');

fprintf('Deleting old figures...\n');
delete(fullfile('figs', 'ip_*.png')); % Delete old figures

if ~exist(OUT_DIR, 'dir') mkdir(OUT_DIR); fprintf('Output directory created: %s\n', OUT_DIR);
else delete(fullfile(OUT_DIR, '*')); fprintf('Output directory already exists. Old files deleted: %s\n', OUT_DIR);
end % Create output directory if it doesn't exist

total_shots = 0;
fprintf('Shots: %s\n', mat2str(shots));
fprintf('\nStarting data retrieval loop...\n');
shot_processing_times = zeros(1, length(shots));
good_shots = [];

for i = 1:length(shots)
    start_time = tic; % Start timer for processing each shot

    shot = shots(i);
    fprintf('\x1b[33mProcessing shot %d [%d/%d]\x1b[0m\n', shot, i, length(shots));

    try 
        %% load liuqe data
        t = liuqe_times(shot); % get times
        t = t(1:DECIMATION:end); % decimate times

        assert(numel(t) > 1, sprintf('Time vector has insufficient elements: t:%s', mat2str(size(t))));
        
        %% run liuqe
        [L, LX, LY] = liuqe(shot,t,'ndom',NDOM,'iterh',ITERH,'iterfrz',ITERFRZ,'iterq',30);

        t = LX.t;
        nt = length(t);
        assert(nt >= MIN_TIME_SAMPLES, 'Not enough time samples');

        % put threshold on ip
        ip = abs(LY.Ip);
        ip_valid = reshape(ip > 1e4, [], 1); % Ip > 10 kA
        assert(sum(ip_valid) >= MIN_TIME_SAMPLES, 'Not enough valid Ip samples');

        % put to zeros points outside of 2 lcfs
        [rry, zzy] = meshgrid(L.ry, L.zy);
        rq0 = reshape(LY.rq(:,end,1:2,:), 32, 2, nt);
        zq0 = reshape(LY.zq(:,end,1:2,:), 32, 2, nt);
        rq1 = reshape(LY.rq(:,end-1,1:2,:), 32, 2, nt);
        zq1 = reshape(LY.zq(:,end-1,1:2,:), 32, 2, nt);
        rq = (0.499 * rq0 + 0.501 * rq1);
        zq = (0.499 * zq0 + 0.501 * zq1);
        dIy = LY.Iy; % original Iy for plotting later
        for ti = 1:nt
            in1 = inpolygon(rry, zzy, rq(:,1,ti), zq(:,1,ti)); % inside lcfs1
            in2 = inpolygon(rry, zzy, rq(:,2,ti), zq(:,2,ti)); % inside lcfs2
            in_domain = in1 | in2;
            LY.Iy(:,:,ti) = LY.Iy(:,:,ti) .* in_domain;
        end
    
        d = ds_prep_save_shot(shot, LX, LY, ip_valid, 1, OUT_DIR);


        total_shots = total_shots + 1;
        good_shots = [good_shots, shot];

        % Quick plot of Fx and Iy at mid time
        close all;
        figure('Position', [0, 1080/4, 1920/2, 1080/2])
        mid = round(numel(d.t)/2); sgtitle(sprintf('Shot %d, t = %.2f s', shot, t(mid)));
        [rr, zz] = meshgrid(L.rx, L.zx);
        [rry, zzy] = meshgrid(L.ry, L.zy);
        subplot(1,3,1); axis equal; hold on;
        Fx = d.Fx(:,:,mid);
        % contourf(rr, zz, Fx, 20, 'LineColor', 'none'); 
        scatter(rr(:), zz(:), 50, Fx(:), 'filled'); colorbar; title('Fx'); 
        % scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled'); 
        plot(rq(:,1,mid), zq(:,1,mid), 'r-', 'LineWidth', 2);
        plot(rq(:,2,mid), zq(:,2,mid), 'g-', 'LineWidth', 2);
        
        subplot(1,3,2); axis equal; hold on; 
        Iy = d.Iy(:,:,mid);
        % contourf(rry, zzy, Iy, 20, 'LineColor', 'none'); 
        scatter(rry(:), zzy(:), 50, Iy(:), 'filled'); colorbar; title('Iy');        
        % scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled');
        plot(rq(:,1,mid), zq(:,1,mid), 'r-', 'LineWidth', 2);
        plot(rq(:,2,mid), zq(:,2,mid), 'g-', 'LineWidth', 2);
        
        subplot(1,3,3); axis equal; hold on; 
        diy = dIy(:,:,mid);
        scatter(rry(:), zzy(:), 50, diy(:), 'filled'); colorbar; title('original Iy');
        % scatter(rq(:,mid), zq(:,mid), 10, 'k', 'filled');
        plot(rq(:,1,mid), zq(:,1,mid), 'r-', 'LineWidth', 2);
        plot(rq(:,2,mid), zq(:,2,mid), 'g-', 'LineWidth', 2);

        saveas(gcf, fullfile('tmp', sprintf('ft_ds_%d.png', shot)));

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

fprintf('\nGood shots processed: %s\n', mat2str(sort(good_shots)));

mdsdisconnect; % Disconnect from MDSplus

fprintf('\nProcessing complete for all shots.\n');
fprintf('Output files saved in: %s\n', OUT_DIR);

% save a random mat in the /NoTivoli/grandin/ directory to show it is finished
save(fullfile('/NoTivoli/grandin/', 'finetune_ds_done.mat'), 'shot_processing_times');