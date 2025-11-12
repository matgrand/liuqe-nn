    
try

    %% current best speeds:
    clear; close all;
    PLOT = false;
    PLOT = true; % set to true to plot results
    n_ctrl_pts = 25; % number of control points
    N = 1000; % number of iterations for inference time test

    % choose the function
    % inference_func = @(phys, r, z) liuqenn_inference(single(phys), single(r), single(z));
    inference_func = @(phys, r, z) liuqenn_inference_mex(single(phys), single(r), single(z));

    % the following dataset can be copied from /home/grandin/repos/liuqe-ml/dss/demo.mat (or it can be generated with prepare_dataset.ipynb)
    try 
        d = load('dss/demo.mat'); % load demo data (should be in DS_DIR/demo.mat , DS_DIR defined in utils.py)
    catch ME
        % create a ds with random data
        fprintf('Could not load demo dataset dss/demo.mat (must be generated with prepare_ds.ipynb),\ngenerating random data instead...\n');
        d = struct(); 
        d.phys = single(randn(100, 134));
        d.pts = single(randn(100, 500, 2));
        d.Fx = single(randn(100, 500));
        d.Iy = single(randn(100, 500));
        d.Br = single(randn(100, 500));
        d.Bz = single(randn(100, 500));
        d.fg = single(randn(100, 4, 65, 28));
    end

    n_examples = size(d.phys, 1);
    max_n_pts = size(d.pts, 2);

    fprintf('Loaded demo data with %d examples and %d points.\n', n_examples, max_n_pts);
    
    NIN = 134;

    %% test with simplified inputs
    rand_i = randi([1, size(d.phys, 1)], 1, 1);
    % rand_i = 1;
    phys = d.phys(rand_i, :); 
    r = d.pts(rand_i, 1:n_ctrl_pts, 1);
    z = d.pts(rand_i, 1:n_ctrl_pts, 2);

    fprintf('phys size: %s, r size: %s, z size: %s\n', ...
        mat2str(size(phys)), mat2str(size(r)), mat2str(size(z)));

    [Fx, Br, Bz] = inference_func(single(phys), single(r), single(z));
    fprintf('Output: Fx %s [%s], Br %s [%s], Bz %s [%s]\n', ...
        mat2str(size(Fx)), mat2str(class(Fx)), ...
        mat2str(size(Br)), mat2str(class(Br)), ...
        mat2str(size(Bz)), mat2str(class(Bz)));


    % print first 5 elements of results
    fprintf('fx_pred  -> [ %s ]\n', num2str(Fx(1:min(8,end)), '%+.4f '));
    fprintf('fx_true  -> [ %s ]\n', num2str(d.Fx(rand_i, 1:min(8,end)), '%+.4f '));
    fprintf('br_pred  -> [ %s ]\n', num2str(Br(1:min(8,end)), '%+.4f '));
    fprintf('br_true  -> [ %s ]\n', num2str(d.Br(rand_i, 1:min(8,end)), '%+.4f '));
    fprintf('bz_pred  -> [ %s ]\n', num2str(Bz(1:min(8,end)), '%+.4f '));
    fprintf('bz_true  -> [ %s ]\n', num2str(d.Bz(rand_i, 1:min(8,end)), '%+.4f '));

    % use timeit to measure inference time
    fprintf('Measuring inference time...\n');
    for n = [1, 10, 20, 25, 30, 100, 128, 256]
        if n >= 100 n_iter = round(N/3);
        else n_iter = N; end

        phys = d.phys(rand_i, :); 
        r = d.pts(rand_i, 1:n, 1);
        z = d.pts(rand_i, 1:n, 2);

        %add noise to the inputs 
        noise_level = single(0.1);
        phys = phys + noise_level * randn(size(phys));
        r = r + noise_level * randn(size(r));
        z = z + noise_level * randn(size(z));

        try
            t1 = timeit(@() inference_func(single(phys), single(r), single(z)));
        catch ME
            t1 = NaN; % in case of error, set time to NaN
        end
        % fprintf('1 sample,      %d control points: %.1f μs\n', n, t * 1e6);
        rand_idxs = randi([1, size(d.phys, 1)], n_iter, 1);
        vphys = d.phys(rand_idxs, :);
        vr = d.pts(rand_idxs, 1:n, 1);
        vz = d.pts(rand_idxs, 1:n, 2);
        try 
            t2 = timeit(@() test_speed(inference_func, vphys, vr, vz, n));
        catch ME
            t2 = NaN; % in case of error, set time to NaN
        end
        % fprintf('%d samples, %d control points: %.1f μs per sample\n', N, n, t * 1e6 / N);
        fprintf('%d pts -> [%.1f μs (1)] [%.1f μs (%d)]\n', ...
            n, t1 * 1e6, t2 * 1e6 / n_iter, n_iter);
    end

    % n_ctrl_pts = 24; % number of control points
    N = 100000;
    fprintf('Testing inference time for %d iterations...', N);
    rand_idxs = randi([1, size(n_examples, 1)], N, 1);
    times = zeros(1, N);
    outs = zeros(N, n_ctrl_pts, 3);
    ttot = tic;
    for i = 1:N
        ri = rand_idxs(i);
        t1 = tic;
        phys_i = single(d.phys(ri, :)+0.1*randn(size(d.phys(ri, :)))); 
        r_i = single(d.pts(ri, 1:n_ctrl_pts, 1)+0.1*randn(size(d.pts(ri, 1:n_ctrl_pts, 1))));
        z_i = single(d.pts(ri, 1:n_ctrl_pts, 2)+0.1*randn(size(d.pts(ri, 1:n_ctrl_pts, 2))));
        [fx, br, bz] = inference_func(phys_i, r_i, z_i);
        times(i) = toc(t1);
    end
    ttot = toc(ttot);
    fprintf(' Done.\n');
    % print times
    fprintf('Inference time -> %.1f ± %.1f [μs] | max %.1f [μs]\n', ...
        mean(times) * 1e6, std(times) * 1e6, max(times) * 1e6);
    fprintf('Total time -> %.1f [s]\n', ttot);
    fprintf('Frequency -> %.1f [μs] | %.1f [Hz]\n', 1e6*ttot/N, N/ttot);
    % assert(false, 'This is a test, please remove this line to continue.');

    %% Iy inference test
    fprintf('Testing Iy inference...\n');
    rand_i = randi([1, size(d.phys, 1)], 1, 1);
    phys = d.phys(rand_i, :); 
    r = d.pts(rand_i, 1:n_ctrl_pts, 1);
    z = d.pts(rand_i, 1:n_ctrl_pts, 2);

    % Iy = liuqe_Iy_nn_inference(single(phys), single(r), single(z));
    Iy = liuqe_Iy_nn_inference_mex(single(phys), single(r), single(z));
    fprintf('Iy size: %s\n', mat2str(size(Iy)));

    fprintf('Output: Iy %s [%s]\n', ...
        mat2str(size(Iy)), mat2str(class(Iy))); 
    fprintf('iy_pred  -> [ %s ]\n', num2str(Iy(1:min(10,end)), '%+4.0f '));
    fprintf('iy_true  -> [ %s ]\n', num2str(d.Iy(rand_i, 1:min(10,end)), '%+4.0f '));

    % test speed for 10_000 iterations
    N = 10000;
    fprintf('Testing Iy inference time for %d iterations...', N);
    rand_idxs = randi([1, size(d.phys, 1)], N, 1);
    times = zeros(1, N);
    ttot = tic;
    for i = 1:N
        ri = rand_idxs(i);
        t1 = tic;
        phys_i = single(d.phys(ri, :)+0.1*randn(size(d.phys(ri, :)))); 
        r_i = single(d.pts(ri, 1:n_ctrl_pts, 1)+0.1*randn(size(d.pts(ri, 1:n_ctrl_pts, 1))));
        z_i = single(d.pts(ri, 1:n_ctrl_pts, 2)+0.1*randn(size(d.pts(ri, 1:n_ctrl_pts, 2))));
        Iy = liuqe_Iy_nn_inference_mex(phys_i, r_i, z_i);
        times(i) = toc(t1);
    end
    ttot = toc(ttot);
    fprintf(' Done.\n');
    % print times
    fprintf('Iy Inference time -> %.1f ± %.1f [μs] | max %.1f [μs]\n', ...
        mean(times) * 1e6, std(times) * 1e6, max(times) * 1e6); 
    fprintf('Total time -> %.1f [s]\n', ttot);
    fprintf('Frequency -> %.1f [μs] | %.1f [Hz]\n', 1e6*ttot/N, N/ttot);

    n_plots = 0;
    rand_idxs =  randi([1, size(d.phys, 1)], 1, n_plots);
    if PLOT
        fprintf('Plotting %d random examples', n_plots);
        for i = 1:n_plots
            rand_i = rand_idxs(i);
            phys = d.phys(rand_i, :); 
            r = d.pts(rand_i, :, 1);
            z = d.pts(rand_i, :, 2);

            [Fx, Br, Bz] = inference_func(single(phys), single(r), single(z));
            Iy = liuqe_Iy_nn_inference(single(phys), single(r), single(z));

            preds = {Fx, Br, Bz, Iy};
            trues = {d.Fx(rand_i, :), d.Br(rand_i, :), d.Bz(rand_i, :), d.Iy(rand_i, :)};
            labels = {'Fx', 'Br', 'Bz', 'Iy'};
            
            fig = figure('Visible', 'off', 'Position', [100, 100, 1000, 1000]);
            % fig = figure('Visible', 'on', 'Position', [100, 100, 1000, 1000]);
            
            for k = 1:4
                ftrue = trues{k};
                fpred = preds{k};
                fmax = max([ftrue(:); fpred(:)]);
                fmin = min([ftrue(:); fpred(:)]);
                diffmap = 100 * abs(ftrue - fpred) ./ (fmax - fmin);
                sz = 20;

                % Plot true value
                subplot(3,4,3*(k-1)+1);
                scatter(r, z, sz, ftrue, 'filled');
                title([labels{k} ' true']);
                xlabel('r'); ylabel('z');
                axis equal; colorbar;
                % clim([fmin fmax]);

                % Plot predicted value
                subplot(3,4,3*(k-1)+2);
                scatter(r, z, sz, fpred, 'filled');
                title([labels{k} ' pred']);
                xlabel('r'); ylabel('z');
                axis equal; colorbar;
                % clim([fmin fmax]);

                % Plot error map
                subplot(3,4,3*(k-1)+3);
                scatter(r, z, sz, diffmap, 'filled');
                title([labels{k} ' Error [%]']);
                xlabel('r'); ylabel('z');
                axis equal; colorbar;
                
                set(gcf, 'Name', labels{k});
            end
            % save figure
            saveas(gcf, sprintf('test/matlab_%s.png', num2str(i)));
            % close figure
            % close(gcf);
            fprintf('.');
        end
        fprintf(' Done.\n');

        % % plot histogram
        % figure('Visible', 'on', 'Position', [100, 100, 800, 600]);
        % histogram(times * 1e6, 100);
        % title('Inference time histogram');
        % xlabel('Time [μs]');
        % ylabel('Count');
        % ylim([0 100]);
        % grid on;
        % % save figure
        % saveas(gcf, 'test/matlab_inference_time.png');
        % close(gcf);
    end

catch ME
    fprintf('Error: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

function t = test_speed(liuqenn_inference, vphys, vr, vz, ncp)
    start = tic;
    n = size(vphys, 1);
    assert(size(vr,1) == n && size(vz,1) == n, ...
        'Input vectors must have the same number of rows.');
    for i = 1:n
        [fx, br, bz] = liuqenn_inference(single(vphys(i,:)), single(vr(i,1:ncp)), single(vz(i,1:ncp)));
    end
    t = toc(start);
end