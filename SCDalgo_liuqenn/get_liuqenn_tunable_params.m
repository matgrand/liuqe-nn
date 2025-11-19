function TP = get_liuqenn_tunable_params()

    % return an empty struct for now (no tunable params)
    TP = struct();

    % general activation flag
    TP.activate = true; 
    
    TP.ton = 0.2; % [s] time after which liuqeml surrogate is activated
    TP.toff = 2.5; % [s] time after which liuqeml surrogate is deactivated
    
end