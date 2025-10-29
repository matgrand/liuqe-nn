function TP = get_liuqeml_tunable_params()

    % return an empty struct for now (no tunable params)
    TP = struct();

    % general activation flag
    TP.activate = true; 
    
    % plasma current threshold activation
    TP.Ip_threshold = 3e4; % [A] % the liuqeml surrogate will be active only when Ip > Ip_threshold

    TP.ton = 0.2; % [s] time after which liuqeml surrogate is activated
    TP.toff = 2.5; % [s] time after which liuqeml surrogate is deactivated
    
end