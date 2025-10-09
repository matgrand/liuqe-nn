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
    
    %% BACKWARD PASS START: Compute dFx/dr and dFx/dz
    
    % The goal is to compute the gradient of Fx w.r.t the inputs r and z.
    % We apply the chain rule backwards from Fx through the network.
    % We assume a scalar loss L = sum(Fx), so the initial gradient dL/dFx is 1.
    
    % 1. Gradient at the output of the network
    % Fx is the first column of head_emb. The gradient for Fx is 1, and 0 for Br, Bz.
    d_head_emb = zeros(size(head_emb), 'like', head_emb);
    d_head_emb(:, 1) = 1.0;

    % 2. Backpropagate through the final activation function of the head network
    d_out2_head = d_head_emb .* d_act(out2_head, net.head_3_beta);

    % 3. Backpropagate through the second linear layer of the head network
    d_act1_head = d_out2_head * net.head_2_weight;

    % 4. Backpropagate through the first activation of the head network
    d_out1_head = d_act1_head .* d_act(out1_head, net.head_1_beta);

    % 5. Backpropagate through the first linear layer of the head network
    d_x = d_out1_head * net.head_0_weight;

    % 6. Backpropagate through the element-wise multiplication
    % x = phys_emb .* pts_embs. We only need the derivative w.r.t pts_embs.
    d_pts_embs = d_x .* phys_emb; % phys_emb is broadcasted

    % 7. Backpropagate through the final activation of the points network
    d_out2 = d_pts_embs .* d_act(out2, net.pts_enc_3_beta);

    % 8. Backpropagate through the second linear layer of the points network
    d_act1 = d_out2 * net.pts_enc_2_weight;

    % 9. Backpropagate through the first activation of the points network
    d_out1 = d_act1 .* d_act(out1, net.pts_enc_1_beta);

    % 10. Backpropagate through the first linear layer of the points network
    d_pts = d_out1 * net.pts_enc_0_weight;

    % 11. Extract the final derivatives w.r.t r and z from d_pts
    % d_pts is an N x 2 matrix where column 1 is dFx/dr and column 2 is dFx/dz.
    % We transpose them to match the shape of the original inputs (1 x N).
    dFx_dr = d_pts(:, 1)';
    dFx_dt = d_pts(:, 2)'; % Note: This is dFx/dz as per the input variable 'z'.
    
    Bz_pinn = dFx_dr ./ (2 * pi * r);  % Bz from derivative
    Br_pinn = -dFx_dt ./ (2 * pi * r); % Br from derivative

    % % Print the PINN-derived magnetic field components
    % fprintf('Br (from derivative):\n');
    % disp(Br_pinn);
    % fprintf('Bz (from derivative):\n');
    % disp(Bz_pinn);
    % fprintf('----------------------------\n');

    % % print Br and Bz from the network output for comparison
    % fprintf('Br (from network output):\n');
    % disp(Br);
    % fprintf('Bz (from network output):\n');
    % disp(Bz);
    % fprintf('----------------------------\n');
    
    % % print mean abs err of Br and Bz from network vs from derivative
    % mean_abs_err_Br = mean(abs(Br - Br_pinn));
    % mean_abs_err_Bz = mean(abs(Bz - Bz_pinn));
    % fprintf('Mean abs error Br: %.6f, Bz: %.6f\n', mean_abs_err_Br, mean_abs_err_Bz);
    Br = Br_pinn;
    Bz = Bz_pinn;

end

%% Activation function: Swish
function y = act(x, beta)
    coder.inline('always');
    % Use bsxfun for element-wise multiplication to ensure
    % compatibility with older MATLAB versions that don't fully
    % support implicit broadcasting in codegen.
    y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
end

%% Derivative of the Swish activation function
function dy_dx = d_act(x, beta)
    coder.inline('always');
    sig_bx = 1 ./ (1 + exp(bsxfun(@times, -beta, x)));
    % Derivative is: sigmoid(beta*x) * (1 + beta * x * (1 - sigmoid(beta*x)))
    dy_dx = sig_bx .* (1 + bsxfun(@times, beta, x) .* (1 - sig_bx));
end