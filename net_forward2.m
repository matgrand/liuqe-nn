function [Fx, Br, Bz] = net_forward2(phys, r, z)
    %#codegen

    % Load the network weights at COMPILE TIME.
    % The 'coder.const' ensures the 'net' struct is embedded as a constant
    % into the final MEX file.
    persistent net;
    if isempty(net)
        % Ensure you are loading a file where ALL fields are single.
        net_const = coder.load('model_weights.mat'); 
        net = coder.const(net_const);
    end

    %% pts net
    pts = [r', z'];
        
    out1 = pts * net.pts_enc_0_weight' + net.pts_enc_0_bias;
    act1 = swish(out1, net.pts_enc_1_beta);
    out2 = act1 * net.pts_enc_2_weight' + net.pts_enc_2_bias;
    pts_embs = swish(out2, net.pts_enc_3_beta);

    %% phys net
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

    out1_phys = phys_norm * net.phys_enc_0_weight' + net.phys_enc_0_bias;
    act1_phys = swish(out1_phys, net.phys_enc_1_beta);
    out2_phys = act1_phys * net.phys_enc_2_weight' + net.phys_enc_2_bias;
    phys_emb = swish(out2_phys, net.phys_enc_3_beta);

    %% head net
    x = phys_emb .* pts_embs; 

    out1_head = x * net.head_0_weight' + net.head_0_bias;
    act1_head = swish(out1_head, net.head_1_beta);
    out2_head = act1_head * net.head_2_weight' + net.head_2_bias;
    head_emb = swish(out2_head, net.head_3_beta);

    Fx = head_emb(:, 1)';
    Br = head_emb(:, 2)';
    Bz = head_emb(:, 3)';
end

% --- Define the Activation Function (Swish) with CORRECTED TYPE ---
function y = swish(x, beta)
    % The literal '1' must be explicitly cast to single to match the type of 'x'.
    y = x ./ (1 + exp(-beta .* x));
end