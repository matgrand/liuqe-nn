function [Fx, Br, Bz] = liuqe_nn_LX(LX, r, z, use_turbo)

    % run liuqenn_inference on standard meq LX struct inputs
    if nargin < 4
        use_turbo = false;
    end
    phys = [LX.Bm; LX.Ff; LX.Ia; LX.Uf; LX.rBt];
    
    Fx = zeros(size(phys, 2), size(r, 1));
    Br = zeros(size(phys, 2), size(r, 1));
    Bz = zeros(size(phys, 2), size(r, 1));
    for i = 1:size(phys, 2)
        phys_i = phys(:, i);
        r_i = r(:);
        z_i = z(:);
        % [Fxi, Bri, Bzi] = liuqenn_inference_mex(single(phys_i)', single(r_i)', single(z_i)');
        [Fxi, Bri, Bzi] = liuqenn_inference(single(phys_i)', single(r_i)', single(z_i)', use_turbo); % we can use the non-mex version here
        Fx(i,:) = Fxi; Br(i,:) = Bri; Bz(i,:) = Bzi;
    end
    phys = phys';
    % figure; grid on; bar(log(abs(phys(1,:)))); title('Bar Plot of phys 1 (LX)');
    % fprintf('phys (LX): %s \n', mat2str(phys(1,:),4));
    assignin('base', 'phys_LX', phys);
end
