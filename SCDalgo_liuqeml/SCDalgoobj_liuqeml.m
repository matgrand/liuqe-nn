function [obj] = SCDalgoobj_liuqeml()

    obj = SCD_algo('SCDalgo_liuqeml');
    
    % Timing of the algorithm
    obj = obj.settiming(-0.5, 1e-4, 3); % same as fmag

    % Fixed parameters init functions
    obj = obj.addfpinitfcn('SCDalgo_liuqeml_loadfp', 'SCDalgo_liuqeml_fp');

    % Tunable parameters
    obj = obj.addtunparamstruct('SCDalgo_liuqeml_tp', @() get_liuqeml_tunable_params());
    obj = obj.addparameter(SCDclass_mdsparscalarint('rtcalgo.liuqeml.params.demo_param', 'DemoParam', 'srcsrv', 'scd', 'srctree', 'scd', 'mdshelp', 'A demo tunable parameter'));
    obj = obj.addparameter(SCDclass_mdsparvectordouble('rtcalgo.liuqeml.params.Ip_threshold','Ip_threshold','srcsrv','scd','srctree','scd','mdshelp','plasma current threshold for liuqeml activation')); 
    obj = obj.addparameter(SCDclass_mdsparvectordouble('rtcalgo.liuqeml.params.toff','toff','srcsrv','scd','srctree','scd','mdshelp','controller switch-off time')); 
    obj = obj.addparameter(SCDclass_mdsparvectordouble('rtcalgo.liuqeml.params.ton' ,'ton' ,'srcsrv','scd','srctree','scd','mdshelp','controller switch-on time'));
    % NOTE: these nodes does not exists on mdsplus yet

    % Buses
    obj = obj.addbus('', @() SCDalgo_liuqeml_loadbus());

    % Tasks

    % Wavegens 

    
end