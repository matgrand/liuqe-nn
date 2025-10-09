% % old version direct calculation of Br and Bz, no backward pass
% function [Fx, Br, Bz] = net_forward2(phys, r, z)
%     %#codegen

%     % Load the network weights at COMPILE TIME.
%     % The 'coder.const' ensures the 'net' struct is embedded as a constant
%     % into the final MEX file.
%     persistent net;
%     if isempty(net)
%         % Ensure you are loading a file where ALL fields are single.
%         % weights_path = 'data/local/weights.mat';
%         % weights_path = 'data/3066835/weights.mat';
% %         weights_path = 'data/3067354/weights.mat'; 
%         weights_path = 'data/3075260/weights.mat';
%         fprintf('Loading weights from %s...', weights_path);
%         net_const = coder.load(weights_path); 
%         net = coder.const(net_const);
%         fprintf('Success\n');
%     end

%     %% pts net
%     pts = [r', z'];
    
%     % start_pts_net = tic;
%     out1 = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
%     act1 = act(out1, net.pts_enc_1_beta);
%     out2 = bsxfun(@plus, act1 * net.pts_enc_2_weight', net.pts_enc_2_bias);
%     pts_embs = act(out2, net.pts_enc_3_beta);
%     % elapsed_pts_net = toc(start_pts_net);
%     % fprintf('Elapsed time for pts net: %.2f micros\n', elapsed_pts_net * 1e6 / size(pts, 1));

%     %% phys net
%     % start_phys_net = tic;   
%     phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

%     out1_phys = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
%     act1_phys = act(out1_phys, net.phys_enc_1_beta);
%     out2_phys = bsxfun(@plus, act1_phys * net.phys_enc_2_weight', net.phys_enc_2_bias);
%     phys_emb = act(out2_phys, net.phys_enc_3_beta);
%     % elapsed_phys_net = toc(start_phys_net);
%     % fprintf('Elapsed time for phys net: %.2f micros\n', elapsed_phys_net * 1e6);

%     %% head net
%     % start_head_net = tic;
%     x = bsxfun(@times, phys_emb, pts_embs);

%     out1_head = bsxfun(@plus, x * net.head_0_weight', net.head_0_bias);
%     act1_head = act(out1_head, net.head_1_beta);
%     out2_head = bsxfun(@plus, act1_head * net.head_2_weight', net.head_2_bias);
%     head_emb = act(out2_head, net.head_3_beta);

%     % elapsed_head_net = toc(start_head_net);
%     % fprintf('Elapsed time for head net: %.2f micros\n', elapsed_head_net * 1e6 / size(pts, 1));

%     Fx = head_emb(:, 1)';
%     Br = head_emb(:, 2)';
%     Bz = head_emb(:, 3)';
% end

% %% Activation function: Swish
% function y = act(x, beta)
%     coder.inline('always');
%     % Use bsxfun for element-wise multiplication to ensure
%     % compatibility with older MATLAB versions that don't fully
%     % support implicit broadcasting in codegen.
%     y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
% end

% % %speed test results:
% % Measuring inference time...
% % 1 pts -> [67.5 μs (1)] [54.5 μs (1000)]
% % 10 pts -> [77.3 μs (1)] [78.3 μs (1000)]
% % 20 pts -> [92.0 μs (1)] [103.7 μs (1000)]
% % 25 pts -> [137.7 μs (1)] [122.3 μs (1000)]
% % 30 pts -> [126.5 μs (1)] [147.2 μs (1000)]
% % 100 pts -> [527.7 μs (1)] [408.5 μs (333)]
% % 128 pts -> [424.4 μs (1)] [439.8 μs (333)]
% % 256 pts -> [737.2 μs (1)] [808.6 μs (333)]
% % Testing inference time for 100000 iterations... Done.
% % Inference time -> 151.2 ± 41.2 [μs] | max 4242.0 [μs]
% % Total time -> 15.3 [s]
% % Frequency -> 152.6 [μs] | 6552.7 [Hz]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new version with optimized forward and backward passes
function [Fx, Br, Bz] = net_forward2(phys, r, z)
    %#codegen

    % Load the network weights at COMPILE TIME.
    % The 'coder.const' ensures the 'net' struct is embedded as a constant
    % into the final MEX file.
    persistent net;
    if isempty(net)
        % Ensure you are loading a file where ALL fields are single.
        % weights_path = 'data/local/weights.mat';
        % weights_path = 'data/3066835/weights.mat';
%         weights_path = 'data/3067354/weights.mat'; 
        weights_path = 'data/3075260/weights.mat';
        fprintf('Loading weights from %s...', weights_path);
        net_const = coder.load(weights_path); 
        net = coder.const(net_const);
        fprintf('Success\n');
    end

    %% --- OPTIMIZED FORWARD PASS ---
    % The forward pass now also stores intermediate sigmoid results 
    % ('sig1', 'sig2', etc.) needed for the backward pass.

    %% pts net
    pts = [r', z'];
    
    out1 = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
    [act1, sig1] = act_with_grad_prep(out1, net.pts_enc_1_beta); % Store sig1
    out2 = bsxfun(@plus, act1 * net.pts_enc_2_weight', net.pts_enc_2_bias);
    [pts_embs, sig2] = act_with_grad_prep(out2, net.pts_enc_3_beta); % Store sig2

    %% phys net
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

    out1_phys = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
    % We don't need the gradients for the phys_net, so we can use the original act function
    act1_phys = act(out1_phys, net.phys_enc_1_beta);
    out2_phys = bsxfun(@plus, act1_phys * net.phys_enc_2_weight', net.phys_enc_2_bias);
    phys_emb = act(out2_phys, net.phys_enc_3_beta);

    %% head net
    x = bsxfun(@times, phys_emb, pts_embs);

    out1_head = bsxfun(@plus, x * net.head_0_weight', net.head_0_bias);
    [act1_head, sig1_head] = act_with_grad_prep(out1_head, net.head_1_beta); % Store sig1_head
    out2_head = bsxfun(@plus, act1_head * net.head_2_weight', net.head_2_bias);
    [head_emb, sig2_head] = act_with_grad_prep(out2_head, net.head_3_beta); % Store sig2_head

    Fx = head_emb(:, 1)';
    % The original Br and Bz are still available if needed, but will be overwritten
    % Br = head_emb(:, 2)'; 
    % Bz = head_emb(:, 3)';
    
    %% --- OPTIMIZED BACKWARD PASS ---
    % This backward pass uses the pre-calculated sigmoid values from the forward pass.
    
    d_head_emb = zeros(size(head_emb), 'like', head_emb);
    d_head_emb(:, 1) = 1.0;

    % Backpropagate using the faster derivative function and stored sigmoids
    d_out2_head = d_head_emb .* d_act_fast(out2_head, head_emb, sig2_head, net.head_3_beta);
    d_act1_head = d_out2_head * net.head_2_weight;
    d_out1_head = d_act1_head .* d_act_fast(out1_head, act1_head, sig1_head, net.head_1_beta);
    d_x = d_out1_head * net.head_0_weight;

    d_pts_embs = bsxfun(@times, d_x, phys_emb);

    d_out2 = d_pts_embs .* d_act_fast(out2, pts_embs, sig2, net.pts_enc_3_beta);
    d_act1 = d_out2 * net.pts_enc_2_weight;
    d_out1 = d_act1 .* d_act_fast(out1, act1, sig1, net.pts_enc_1_beta);
    d_pts = d_out1 * net.pts_enc_0_weight;

    dFx_dr = d_pts(:, 1)';
    dFx_dt = d_pts(:, 2)';
    
    % Calculate final values from derivatives
    Bz_pinn = dFx_dr ./ (2 * pi * r);
    Br_pinn = -dFx_dt ./ (2 * pi * r);

    Br = Br_pinn;
    Bz = Bz_pinn;

end

%% --- OPTIMIZED ACTIVATION FUNCTIONS ---

% Original activation function for paths where no gradient is needed (phys_net)
function y = act(x, beta)
    coder.inline('always');
    y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
end

% New function: Computes activation AND prepares the sigmoid for the backward pass
function [y, sig_bx] = act_with_grad_prep(x, beta)
    coder.inline('always');
    sig_bx = 1 ./ (1 + exp(bsxfun(@times, -beta, x)));
    y = x .* sig_bx;
end

% New function: Computes derivative much faster using the pre-computed sigmoid
function dy_dx = d_act_fast(x, y, sig_bx, beta)
    coder.inline('always');
    % y is the activation result: x .* sig_bx
    % The full derivative is: sig_bx + beta * y * (1 - sig_bx)
    dy_dx = sig_bx + bsxfun(@times, beta, y) .* (1 - sig_bx);
end

% % %speed test results:
% Measuring inference time...
% 1 pts -> [80.5 μs (1)] [61.4 μs (1000)]
% 10 pts -> [84.5 μs (1)] [92.4 μs (1000)]
% 20 pts -> [120.2 μs (1)] [168.1 μs (1000)]
% 25 pts -> [141.6 μs (1)] [150.2 μs (1000)]
% 30 pts -> [156.0 μs (1)] [164.6 μs (1000)]
% 100 pts -> [463.0 μs (1)] [460.3 μs (333)]
% 128 pts -> [802.7 μs (1)] [550.8 μs (333)]
% 256 pts -> [991.5 μs (1)] [984.6 μs (333)]
% Testing inference time for 100000 iterations... Done.
% Inference time -> 182.7 ± 42.2 [μs] | max 2642.0 [μs]
% Total time -> 18.4 [s]
% Frequency -> 184.2 [μs] | 5428.4 [Hz]