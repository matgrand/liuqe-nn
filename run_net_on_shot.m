function run_net_on_shot(shot_number, save_dir)
    fprintf('Running network on shot %d, saving in %s\n', shot_number, save_dir);
    try
        % addpath([pwd '/onnx_net_forward']);
        % model_path = [pwd '/onnx_net_forward/net.onnx'];

        % % first call to load the model
        % liuqenn_inference_mex(model_path);


        % load shot
        shot_file_path = fullfile('test_shots', sprintf('%d_cache.mat', shot_number));

        d = load(shot_file_path);
        fprintf('Loaded shot data from: %s\n', shot_file_path);

        phys = pack_inputs(d.Bm, d.Ff, d.Ft, d.Ia, d.Ip, d.Iu, d.rBt);
        fprintf('Packed inputs size: %s\n', mat2str(size(phys)));

        nt = size(phys,2); % number of time points

        %load tcv grid
        g = load('tcv_params/grid.mat');
        gr = g.r;
        gz = g.z;

        % create a meshgrid from r and z
        [r, z] = meshgrid(gr, gz);
        %flatten r and z
        r = single(r(:));
        z = single(z(:));

        % control points
        rc = [0.7038, 0.6722, 0.6516, 0.7108, 0.9376, 1.0843, 1.0931, 0.9414, 0.8023, 0.6240];
        zc = [-0.1195, 0.1285, 0.3775, 0.6159, 0.6127, 0.4150, 0.1691, -0.0246, -0.7500, -0.1229];
        rq = rc;
        zq = zc;
        n_ctrl_pts = length(rc);
        assert(n_ctrl_pts == length(zc), 'Control points rc and zc must have the same length');

        Fx = zeros(65*28, nt); % preallocate Fx
        Br = zeros(65*28, nt); % preallocate Br
        Bz = zeros(65*28, nt); % preallocate Bz
        
        for i = 1:nt % loop over time points
            % [Fx(:,i), Br(:,i), Bz(:,i)] = liuqenn_inference_mex(phys(:, i), r, z);
            [Fx(:,i), Br(:,i), Bz(:,i)] = liuqenn_inference(phys(:, i)', r', z');
        end

        Fxc = zeros(n_ctrl_pts, nt); % preallocate Fx at control points
        Brc = zeros(n_ctrl_pts, nt); % preallocate Br at control points
        Bzc = zeros(n_ctrl_pts, nt); % preallocate Bz at control points
        for i = 1:nt % loop over time points
            % [Fxc(:,i), Brc(:,i), Bzc(:,i)] = liuqenn_inference_mex(phys(:, i), rc', zc');
            [Fxc(:,i), Brc(:,i), Bzc(:,i)] = liuqenn_inference(phys(:, i)', rc, zc);
        end

        Fxq = zeros(n_ctrl_pts, nt); % preallocate Fx at control points
        Brq = zeros(n_ctrl_pts, nt); % preallocate Br at control points
        Bzq = zeros(n_ctrl_pts, nt); % preallocate Bz at control points
        for i = 1:nt % loop over time points
            % [Fxq(:,i), Brq(:,i), Bzq(:,i)] = liuqenn_inference_mex(phys(:, i), rq', zq');
            [Fxq(:,i), Brq(:,i), Bzq(:,i)] = liuqenn_inference(phys(:, i)', rq, zq);
        end

        % reshape results to match the grid size
        Fx = reshape(Fx, 65, 28, nt);
        Br = reshape(Br, 65, 28, nt);
        Bz = reshape(Bz, 65, 28, nt);

        assert(all(size(Fx) == size(d.Fx)), 'Fx has wrong size');
        assert(all(size(Br) == size(d.Br)), 'Br has wrong size');
        assert(all(size(Bz) == size(d.Bz)), 'Bz has wrong size');

        % save results in a .mat file
        save_file = fullfile(save_dir, sprintf('%d_net.mat', shot_number));
        save(save_file, 'Fx', 'Br', 'Bz', 'Fxc', 'Brc', 'Bzc', 'Fxq', 'Brq', 'Bzq', 'rq', 'zq');
        fprintf('Saved results to: %s\n', save_file);

    catch ME
        fprintf('Error: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end

    function inputs = pack_inputs(Bm, Ff, Ft, Ia, Ip, Iu, rBt)
        % Pack inputs into a single array 
        % NOTE: the order is very important here, it should match the order in which the network was trained, see utils.py -> INPUT_NAMES
        inputs = [Bm; Ff; Ft; Ia; Ip; Iu; rBt];
        % Ensure inputs are in single precision
        inputs = single(inputs);
    end
end
