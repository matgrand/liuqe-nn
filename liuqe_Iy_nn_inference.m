% old version direct calculation of Br and Bz, no backward pass
function [Iy] = liuqe_Iy_nn_inference(phys, r, z)
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
        % weights_path = 'pretrained/weights_turbosmall_l1.mat'; 

        fprintf('Loading weights from %s...   (make sure to run the script pretrained/json2mat.m to convert JSON weights to MAT files)\n', weights_path);
        net_const = coder.load(weights_path); 
        net = coder.const(net_const);
        
        fields = fieldnames(net);
        for i = 1:numel(fields)
            sz = size(net.(fields{i}));
            sz_str = sprintf('%dx', sz);
            sz_str = sz_str(1:end-1); % remove trailing 'x'
            fprintf('net.%s: [%s]\n', fields{i}, sz_str);
        end
        
        fprintf('Success\n');
    end

    %% pts 
    pts = [r', z'];
    
    pts_embs = bsxfun(@plus, pts * net.pts_enc_0_weight', net.pts_enc_0_bias);
    pts_embs = act(pts_embs, net.pts_enc_1_beta);
    pts_embs = bsxfun(@plus, pts_embs * net.pts_enc_2_weight', net.pts_enc_2_bias);
    pts_embs = act(pts_embs, net.pts_enc_3_beta);

    %% phys net
    phys_norm = (phys - net.x_mean_std(1,:)) ./ net.x_mean_std(2,:);

    phys_emb = bsxfun(@plus, phys_norm * net.phys_enc_0_weight', net.phys_enc_0_bias);
    phys_emb = act(phys_emb, net.phys_enc_1_beta);
    phys_emb = bsxfun(@plus, phys_emb * net.phys_enc_2_weight', net.phys_enc_2_bias);
    phys_emb = act(phys_emb, net.phys_enc_3_beta);

    %% position + physics aggregation
    x = bsxfun(@times, phys_emb, pts_embs);

    %% Iy net
    iy1 = bsxfun(@plus, x * net.iy_head_0_weight', net.iy_head_0_bias);
    aiy1 = act(iy1, net.iy_head_1_beta);
    iy2 = bsxfun(@plus, aiy1 * net.iy_head_2_weight', net.iy_head_2_bias);
    Iy = act(iy2, net.iy_head_3_beta)';
    Iy = Iy * 0.00046884178333783666 / (4*pi*1e-7); % scale by grid_area/mu0

end

function y = act(x, beta)
    coder.inline('always');
    y = x ./ (1 + exp(bsxfun(@times, -beta, x)));
end