function [Fx, Br, Bz] = liuqenn_inference_select(phys, r, z, weights)
    %#codegen
    [Fx, Br, Bz] = liuqenn_inference(phys, r, z, weights);
end
