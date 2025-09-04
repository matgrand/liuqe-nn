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

    %% pts net
    pts = [r', z'];
    
    % start_pts_net = tic;
    out1 = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
    act1 = act(out1, net.pts_enc_1_beta);
    out2 = bsxfun(@plus, act1 * net.pts_enc_2_weight', net.pts_enc_2_bias);
    pts_embs = act(out2, net.pts_enc_3_beta);
    % elapsed_pts_net = toc(start_pts_net);
    % fprintf('Elapsed time for pts net: %.2f micros\n', elapsed_pts_net * 1e6 / size(pts, 1));

    %% phys net
    % start_phys_net = tic;   
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

    out1_phys = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
    act1_phys = act(out1_phys, net.phys_enc_1_beta);
    out2_phys = bsxfun(@plus, act1_phys * net.phys_enc_2_weight', net.phys_enc_2_bias);
    phys_emb = act(out2_phys, net.phys_enc_3_beta);
    % elapsed_phys_net = toc(start_phys_net);
    % fprintf('Elapsed time for phys net: %.2f micros\n', elapsed_phys_net * 1e6);

    %% head net
    % start_head_net = tic;
    x = bsxfun(@times, phys_emb, pts_embs);

    out1_head = bsxfun(@plus, x * net.head_0_weight', net.head_0_bias);
    act1_head = act(out1_head, net.head_1_beta);
    out2_head = bsxfun(@plus, act1_head * net.head_2_weight', net.head_2_bias);
    head_emb = act(out2_head, net.head_3_beta);

    % elapsed_head_net = toc(start_head_net);
    % fprintf('Elapsed time for head net: %.2f micros\n', elapsed_head_net * 1e6 / size(pts, 1));

    Fx = head_emb(:, 1)';
    Br = head_emb(:, 2)';
    Bz = head_emb(:, 3)';
end

%% Activation function: Swish
function y = act(x, beta)
    coder.inline('always');
    % Use bsxfun for element-wise multiplication to ensure
    % compatibility with older MATLAB versions that don't fully
    % support implicit broadcasting in codegen.
    y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
end
