function [outBus, memBus] = liuqeml(inBus, SCDalgo_liuqeml_fp,SCDalgo_liuqeml_tp)
    % this is a copy of the function inside the SCDalgo_liuqeml.slx simulink block
    
    % inBus: [Bm, Ff, Ia, Ip, Ft, Iu, rBt, rn, zn]
    % outBus: [Fn, Brn, Bzn]
    % memBus: [*inBus, *tpBus] (just log everything)

    % addtional inputs: 
    %  - fixed params = SCDalgo_liuqeml_fp
    %  - tunable params = SCDalgo_liuqeml_tp

    FP = SCDalgo_liuqeml_fp;
    TP = SCDalgo_liuqeml_tp;
    
    % do something with tunable params (like update net weights)
    % (for now there are no tunable params)
        
    % unpack inputs (that I would like)
    Bm = inBus.Bm; Ff = inBus.Ff; Ia = inBus.Ia; Ip = inBus.Ip;
    Ft = inBus.Ft; Iu = inBus.Iu; rBt = inBus.rBt;
    rn = inBus.rn; zn = inBus.zn;
    
    % % unpack the probably actual inputs
    % Bm = inBus.Bm; % should be in liuin
    % Ff = inBus.Ff; % should be in liuin
    % rBt = inBus.rBt; % should be in liuin
    % Ia = inBus.Ipf; % ? <IPF> from liuin, may require remapping ->  Ia = [LY.Ia(cpar.iEL) LY.Ia(cpar.iFL) LY.Ia(cpar.iGL) LY.Ia(cpar.iOHL)];
    % Ip = inBus.Ip; % ?  
    % Ft = inBus.Ft; % ?
    % Iu = inBus.Iu; % manipulations from Uf ->  Iu = -LY.Uf./cpar.Rs;
    % rn = inBus.rn; zn = inBus.zn; % points

    % run network inference
    phys = [Bm; Ff; Ft; Ia; Ip; Iu; rBt];
    [Fxn,Brn,Bzn] = net_forward2_mex(phys,rn,zn);

    % assign outputs
    outBus = struct();
    outBus.Fn = Fxn; outBus.Brn = Brn; outBus.Bzn = Bzn;

    % log everything in memBus
    memBus = struct();
    memBus.Fn = Fxn; memBus.Brn = Brn; memBus.Bzn = Bzn;
    memBus.Bm = Bm; memBus.Ff = Ff; memBus.Ia = Ia; memBus.Ip = Ip;
    memBus.Ft = Ft; memBus.Iu = Iu; memBus.rBt = rBt;
    memBus.rn = rn; memBus.zn = zn;

end