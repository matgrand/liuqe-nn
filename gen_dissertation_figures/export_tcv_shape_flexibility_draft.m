% Export the twelve LIUQE equilibria shown in tcv_shape_flexibility_draft.png.
% Run on the TCV server from any directory, then copy data/tcv_shape_flexibility_draft/*.mat to the same repository path locally.
shots = [9849, 6010, 5650, 8856, 11962, 11368, 18548, 6442, 10159, 19373, 11928, 36151, 87611];
times = [0.5, 0.8, 0.8, 0.29, 1.0, 0.65, 1.5, 0.5, 0.012, 0.42, 0.67, 0.457, 0.8];
repo_dir = fileparts(fileparts(mfilename('fullpath')));
output_dir = fullfile(repo_dir, 'data', 'tcv_shape_flexibility_draft');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

try
    mdsconnect('tcvdata.epfl.ch');
catch err
    warning('MDSplus connection: %s', err.message);
end

for k = 1:numel(shots)
    shot = shots(k); requested_t = times(k); loaded = false;
    filename = sprintf('shot_%d_t_%0.3f.mat', shot, requested_t);
    output_path = fullfile(output_dir, filename);
    try
        [L, LY] = mds2meq(shot, 'LIUQE.M');
        if ~isfield(LY, 't') || isempty(LY.t), error('Stored equilibrium has no time vector.'); end
        [dt, idx] = min(abs(LY.t(:) - requested_t));
        if ~isfinite(dt) || dt > 0.0015, error('Nearest stored equilibrium is %.4f s from the requested time.', dt); end
        [r, z, psi, psi_axis, psi_boundary] = shape_slice(L, LY, idx);
        actual_t = LY.t(idx); source = 'LIUQE.M'; loaded = true;
    catch err
        fprintf(2, 'Stored equilibrium #%d at %.3f s unavailable: %s\n', shot, requested_t, err.message);
    end
    if ~loaded
        try
            [L, LX, LY] = liuqe(shot, requested_t);
            [r, z, psi, psi_axis, psi_boundary] = shape_slice(L, LY, 1);
            actual_t = requested_t;
            if isfield(LX, 't') && ~isempty(LX.t), actual_t = LX.t(1); end
            source = 'liuqe'; loaded = true;
        catch err
            fprintf(2, 'FAILED #%d at %.3f s: %s\n', shot, requested_t, err.message);
        end
    end
    if ~loaded, continue; end
    try
        schema_version = 1;
        save(output_path, 'schema_version', 'shot', 'requested_t', 'actual_t', 'source', 'r', 'z', 'psi', 'psi_axis', 'psi_boundary');
        fprintf('Saved %s (%s, actual t=%.4f s)\n', output_path, source, actual_t);
    catch err
        fprintf(2, 'Could not save #%d at %.3f s: %s\n', shot, requested_t, err.message);
    end
end

function [r, z, psi, psi_axis, psi_boundary] = shape_slice(L, LY, idx)
    r = L.rx(:); z = L.zx(:);
    if ~isfield(LY, 'Fx') || ~isfield(LY, 'FA') || ~isfield(LY, 'FB')
        error('Equilibrium is missing Fx, FA, or FB.');
    end
    psi = LY.Fx(:, :, idx);
    if isequal(size(psi), [numel(r), numel(z)]), psi = psi.'; end
    if ~isequal(size(psi), [numel(z), numel(r)]), error('Unexpected flux-map dimensions.'); end
    psi_axis = LY.FA(1, idx); psi_boundary = LY.FB(1, idx);
    if any(~isfinite(psi(:))) || ~isfinite(psi_axis) || ~isfinite(psi_boundary) || abs(psi_axis - psi_boundary) < eps
        error('Equilibrium has invalid flux values.');
    end
end
