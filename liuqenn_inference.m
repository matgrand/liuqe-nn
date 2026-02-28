% old version direct calculation of Br and Bz, no backward pass
function [Fx, Br, Bz] = liuqenn_inference(phys, r, z)
    %#codegen

    % Load the network weights at COMPILE TIME.
    % The 'coder.const' ensures the 'net' struct is embedded as a constant
    % into the final MEX file.
    persistent net;
    if isempty(net)
        %% NOTE: to generate the MAT file from JSON weights, run the script: pretrained/json2mat.m
        % uncomment one of the following weight files to use different models
        weights_path = 'pretrained/weights_small_l1.mat'; % <- best
        % weights_path = 'pretrained/weights_small_l2.mat'; 
        % weights_path = 'pretrained/weights_turbo_l1.mat'; 

        fprintf('Loading weights from %s...   (make sure to run the script pretrained/json2mat.m to convert JSON weights to MAT files)\n', weights_path);
        net_const = coder.load(weights_path); 
        net = coder.const(net_const);
        
        % fields = fieldnames(net);
        % for i = 1:numel(fields)
        %     sz = size(net.(fields{i}));
        %     sz_str = sprintf('%dx', sz);
        %     sz_str = sz_str(1:end-1); % remove trailing 'x'
        %     fprintf('net.%s: [%s]\n', fields{i}, sz_str);
        % end
        
        % fprintf('Success\n');
    end

    % make sure all r > 0
    r(r <= 0) = 1e-3;

    % [Fx, Br, Bz] = direct_net_forward(phys, r, z, net); % standard version
    [Fx, Br, Bz] = grad_net_forward(phys, r, z, net); % gradient version
end

function [Fx, Br, Bz] = direct_net_forward(phys, r, z, net)
    %#codegen
    %% pts 
    pts = [r', z'];
    
    pts_embs = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
    pts_embs = act(pts_embs, net.pts_enc_1_beta);
    pts_embs = bsxfun(@plus, pts_embs * net.pts_enc_2_weight', net.pts_enc_2_bias);
    pts_embs = act(pts_embs, net.pts_enc_3_beta);
    % pts_embs = bsxfun(@plus, pts_embs * net.pts_enc_4_weight', net.pts_enc_4_bias);
    % pts_embs = act(pts_embs, net.pts_enc_5_beta);

    %% phys net
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

    phys_emb = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
    phys_emb = act(phys_emb, net.phys_enc_1_beta);
    phys_emb = bsxfun(@plus, phys_emb * net.phys_enc_2_weight', net.phys_enc_2_bias);
    phys_emb = act(phys_emb, net.phys_enc_3_beta);
    % phys_emb = bsxfun(@plus, phys_emb * net.phys_enc_4_weight', net.phys_enc_4_bias);
    % phys_emb = act(phys_emb, net.phys_enc_5_beta);

    %% position + physics aggregation
    x = bsxfun(@times, phys_emb, pts_embs);
    
    %% Fx net
    f = bsxfun(@plus, x * net.fx_head_0_weight', net.fx_head_0_bias);
    f = act(f, net.fx_head_1_beta);
    f = bsxfun(@plus, f * net.fx_head_2_weight', net.fx_head_2_bias);
    f = act(f, net.fx_head_3_beta);
    % f = bsxfun(@plus, f * net.fx_head_4_weight', net.fx_head_4_bias);
    % f = act(f, net.fx_head_5_beta);

    % %% Iy net
    % fx1 = bsxfun(@plus, x * net.iy_head_0_weight', net.iy_head_0_bias);
    % afx1 = act(fx1, net.iy_head_1_beta);
    % fx2 = bsxfun(@plus, afx1 * net.iy_head_2_weight', net.iy_head_2_bias);
    % Iy = act(fx2, net.iy_head_3_beta);
    
    %% Br Bz net
    b = bsxfun(@plus, x * net.brz_head_0_weight', net.brz_head_0_bias);
    b = act(b, net.brz_head_1_beta);
    b = bsxfun(@plus, b * net.brz_head_2_weight', net.brz_head_2_bias);
    b = act(b, net.brz_head_3_beta);
    % b = bsxfun(@plus, b * net.brz_head_4_weight', net.brz_head_4_bias);
    % b = act(b, net.brz_head_5_beta);
    
    Fx = f'; % assign Fx
    Br = b(:, 1)'; % assign Br
    Bz = b(:, 2)'; % assign Bz
end

function [Fx, Br, Bz] = grad_net_forward(phys, r, z, net)

    %% pts net
    pts = [r', z'];
    
    p0 = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
    [ap1, sp1] = act_with_grad_prep(p0, net.pts_enc_1_beta); % Store sp1
    p2 = bsxfun(@plus, ap1 * net.pts_enc_2_weight', net.pts_enc_2_bias);
    [pts_embs, sp2] = act_with_grad_prep(p2, net.pts_enc_3_beta); % Store sp2

    %% phys net (we dont need grad wrt phys)
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);
    phys_emb = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
    phys_emb = act(phys_emb, net.phys_enc_1_beta);
    phys_emb = bsxfun(@plus, phys_emb * net.phys_enc_2_weight', net.phys_enc_2_bias);
    phys_emb = act(phys_emb, net.phys_enc_3_beta);
    % phys_emb = bsxfun(@plus, phys_emb * net.phys_enc_4_weight', net.phys_enc_4_bias);
    % phys_emb = act(phys_emb, net.phys_enc_5_beta);

    %% position + physics aggregation
    x = bsxfun(@times, phys_emb, pts_embs);
    
    %% Fx net
    fx1 = bsxfun(@plus, x * net.fx_head_0_weight', net.fx_head_0_bias);
    [afx1, sfx1] = act_with_grad_prep(fx1, net.fx_head_1_beta); % Store sfx1
    fx2 = bsxfun(@plus, afx1 * net.fx_head_2_weight', net.fx_head_2_bias);
    [fx3, sfx3] = act_with_grad_prep(fx2, net.fx_head_3_beta); % Store sfx3

    Fx = fx3'; % assign Fx

    %% backward pass to get dFx/dr and dFx/dz
    dfx2 = 1.0 .* d_act_fast(fx2, fx3, sfx3, net.fx_head_3_beta);
    dafx1 = dfx2 * net.fx_head_2_weight;
    dfx1 = dafx1 .* d_act_fast(fx1, afx1, sfx1, net.fx_head_1_beta);
    dx = dfx1 * net.fx_head_0_weight;

    d_pts_embs = bsxfun(@times, dx, phys_emb);

    dp2 = d_pts_embs .* d_act_fast(p2, pts_embs, sp2, net.pts_enc_3_beta);
    dap1 = dp2 * net.pts_enc_2_weight;
    dp1 = dap1 .* d_act_fast(p0, ap1, sp1, net.pts_enc_1_beta);
    d_pts = dp1 * net.pts_enc_0_weight;

    dFx_dr = d_pts(:, 1)';
    dFx_dt = d_pts(:, 2)';
    
    % Calculate final values from derivatives
    Bz = dFx_dr ./ (2 * pi * r);
    Br = -dFx_dt ./ (2 * pi * r);
end

%% activation functions

% Original activation function for paths where no gradient is needed (phys_net)
function y = act(x, beta)
    coder.inline('always');
    y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
end

% activation function: Computes activation AND prepares the sigmoid for the backward pass
function [y, sig_bx] = act_with_grad_prep(x, beta)
    coder.inline('always');
    sig_bx = 1 ./ (1 + exp(bsxfun(@times, -beta, x)));
    y = x .* sig_bx;
end

% activation function: Computes derivative using the pre-computed sigmoid
function dy_dx = d_act_fast(x, y, sig_bx, beta)
    coder.inline('always');
    % y is the activation result: x .* sig_bx
    % The full derivative is: sig_bx + beta * y * (1 - sig_bx)
    dy_dx = sig_bx + bsxfun(@times, beta, y) .* (1 - sig_bx);
end
