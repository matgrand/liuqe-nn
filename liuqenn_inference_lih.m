function [Fx, dFx_dr, dFx_dt, d2Fx_dr2, d2Fx_dz2] = liuqenn_inference_lih(phys, r, z)
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

    [Fx, dFx_dr, dFx_dt, d2Fx_dr2, d2Fx_dz2] = grad_net_forward(phys, r, z, net); % gradient version
end

function [Fx, dFx_dr, dFx_dt, d2Fx_dr2, d2Fx_dz2] = grad_net_forward(phys, r, z, net)

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

    %% Forward-mode AD: first and second derivatives w.r.t. r and z
    % For a linear layer p = t*W' + b: dp/dr = W(:,1)', d²p/dr² = 0
    % For a swish layer y = swish(x,β): dy = d1*dx, d²y = d2*dx² + d1*d²x
    %   where d1 = d_act_fast, d2 = d2_act (see helpers below)

    % Layer pts_enc_0 (linear)
    dr_p0 = net.pts_enc_0_weight(:, 1)';  % [1×h0], broadcast over n
    dz_p0 = net.pts_enc_0_weight(:, 2)';  % [1×h0]

    % Layer pts_enc_1 (swish)
    d1_a1  = d_act_fast(p0, ap1, sp1, net.pts_enc_1_beta);  % [n×h0]
    d2_a1  = d2_act(p0, sp1, net.pts_enc_1_beta);           % [n×h0]
    dr_ap1  = d1_a1 .* dr_p0;
    dz_ap1  = d1_a1 .* dz_p0;
    dr2_ap1 = d2_a1 .* (dr_p0 .^ 2);   % d²p0/dr²=0 term vanishes
    dz2_ap1 = d2_a1 .* (dz_p0 .^ 2);

    % Layer pts_enc_2 (linear)
    dr_p2  = dr_ap1  * net.pts_enc_2_weight';
    dz_p2  = dz_ap1  * net.pts_enc_2_weight';
    dr2_p2 = dr2_ap1 * net.pts_enc_2_weight';
    dz2_p2 = dz2_ap1 * net.pts_enc_2_weight';

    % Layer pts_enc_3 (swish)
    d1_pe  = d_act_fast(p2, pts_embs, sp2, net.pts_enc_3_beta);  % [n×h2]
    d2_pe  = d2_act(p2, sp2, net.pts_enc_3_beta);                % [n×h2]
    dr_pe  = d1_pe .* dr_p2;
    dz_pe  = d1_pe .* dz_p2;
    dr2_pe = d2_pe .* (dr_p2 .^ 2) + d1_pe .* dr2_p2;
    dz2_pe = d2_pe .* (dz_p2 .^ 2) + d1_pe .* dz2_p2;

    % Aggregation: x = phys_emb .* pts_embs  (phys_emb constant w.r.t. r, z)
    dr_x  = phys_emb .* dr_pe;
    dz_x  = phys_emb .* dz_pe;
    dr2_x = phys_emb .* dr2_pe;
    dz2_x = phys_emb .* dz2_pe;

    % Layer fx_head_0 (linear)
    dr_fx1  = dr_x  * net.fx_head_0_weight';
    dz_fx1  = dz_x  * net.fx_head_0_weight';
    dr2_fx1 = dr2_x * net.fx_head_0_weight';
    dz2_fx1 = dz2_x * net.fx_head_0_weight';

    % Layer fx_head_1 (swish)
    d1_a2  = d_act_fast(fx1, afx1, sfx1, net.fx_head_1_beta);  % [n×hfx]
    d2_a2  = d2_act(fx1, sfx1, net.fx_head_1_beta);            % [n×hfx]
    dr_afx1  = d1_a2 .* dr_fx1;
    dz_afx1  = d1_a2 .* dz_fx1;
    dr2_afx1 = d2_a2 .* (dr_fx1 .^ 2) + d1_a2 .* dr2_fx1;
    dz2_afx1 = d2_a2 .* (dz_fx1 .^ 2) + d1_a2 .* dz2_fx1;

    % Layer fx_head_2 (linear)
    dr_fx2  = dr_afx1  * net.fx_head_2_weight';
    dz_fx2  = dz_afx1  * net.fx_head_2_weight';
    dr2_fx2 = dr2_afx1 * net.fx_head_2_weight';
    dz2_fx2 = dz2_afx1 * net.fx_head_2_weight';

    % Layer fx_head_3 (swish)
    d1_fx3v = d_act_fast(fx2, fx3, sfx3, net.fx_head_3_beta);  % [n×1]
    d2_fx3v = d2_act(fx2, sfx3, net.fx_head_3_beta);           % [n×1]
    dr_fx3  = d1_fx3v .* dr_fx2;
    dz_fx3  = d1_fx3v .* dz_fx2;
    dr2_fx3 = d2_fx3v .* (dr_fx2 .^ 2) + d1_fx3v .* dr2_fx2;
    dz2_fx3 = d2_fx3v .* (dz_fx2 .^ 2) + d1_fx3v .* dz2_fx2;

    dFx_dr   = dr_fx3';   % [1×n]
    dFx_dt   = dz_fx3';   % [1×n]
    d2Fx_dr2 = dr2_fx3';  % [1×n]
    d2Fx_dz2 = dz2_fx3';  % [1×n]

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

% Second derivative of swish: d²(x·σ(βx))/dx² = β·σ·(1-σ)·(2 + β·x·(1-2σ))
function d2y_dx2 = d2_act(x, sig_bx, beta)
    coder.inline('always');
    d2y_dx2 = bsxfun(@times, beta, sig_bx) .* (1 - sig_bx) .* ...
              (2 + bsxfun(@times, beta, x) .* (1 - 2*sig_bx));
end
