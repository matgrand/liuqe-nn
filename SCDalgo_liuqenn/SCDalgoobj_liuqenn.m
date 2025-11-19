function [obj] = SCDalgoobj_liuqenn()

    obj = SCD_algo('SCDalgo_liuqenn');

    % Timing of the algorithm
    obj = obj.settiming(-0.5, 1e-4, 3); % same as fmag

    % Fixed parameters init functions
    obj = obj.addfpinitfcn('SCDalgo_liuqenn_loadfp', 'SCDalgo_liuqenn_fp');

    % Tunable parameters
    obj = obj.addtunparamstruct('SCDalgo_liuqenn_tp', @() get_liuqenn_tunable_params());
    obj = obj.addparameter(SCDclass_mdsparvectordouble('rtcalgo.liuqenn.params.toff','toff','srcsrv','scd','srctree','scd','mdshelp','controller switch-off time')); 
    obj = obj.addparameter(SCDclass_mdsparvectordouble('rtcalgo.liuqenn.params.ton' ,'ton' ,'srcsrv','scd','srctree','scd','mdshelp','controller switch-on time'));
    obj = obj.addparameter(SCDclass_mdsparvectorbool('rtcalgo.liuqenn.params.activate','activate','srcsrv','scd','srctree','scd','mdshelp','Flag to activate liuqenn'));


    % NOTE: these nodes does not exists on mdsplus yet

    % Buses
    obj = obj.addbus('', @() SCDalgo_liuqenn_loadbus());

    % Tasks

    % Wavegens 
    obj=obj.addwavegen(SCDclass_mdswgsigarray1('rtcalgo.liuqenn.waves.rn_%03d','sc_ri',1:25,'srcsrv','scd','srctree','scd'));
    obj=obj.addwavegen(SCDclass_mdswgsigarray1('rtcalgo.liuqenn.waves.zn_%03d','sc_zi',1:25,'srcsrv','scd','srctree','scd'));
    
end