function [name, bus] = SCDalgo_liuqenn_loadbus()

    %% Fixed parameters
    FP = SCDalgo_liuqenn_loadfp();

    %% inBus t + [BM, FF, IA, UF, RBT] + [rn, zn]
    ii = 1;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 't';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [1 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = '';

    ii = 2;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Bm';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.nBm) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = '';

    ii = 3;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Ff';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.nFf) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 4;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Ia';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.nIa) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 5;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Uf';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.nUf) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 6;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'rBt';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.nrBt) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 7;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'rn';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.npts) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];     
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 8;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'zn';
    tmp{1}(ii, 1).Complexity = 'real';  
    tmp{1}(ii, 1).Dimensions = [double(FP.npts) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = '';     

    inBus = Simulink.Bus;
    inBus.Description = '';
    inBus.DataScope = 'Auto';
    inBus.HeaderFile = '';
    inBus.Alignment = -1;
    inBus.Elements = tmp{1};    
    clear tmp;

    %% outBus
    ii = 1;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Fxn';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.npts) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 2;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Brn';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.npts) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';    
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    ii = 3;
    tmp{1}(ii, 1) = Simulink.BusElement;
    tmp{1}(ii, 1).Name = 'Bzn';
    tmp{1}(ii, 1).Complexity = 'real';
    tmp{1}(ii, 1).Dimensions = [double(FP.npts) 1];
    tmp{1}(ii, 1).DataType = 'single';
    tmp{1}(ii, 1).Min = [];
    tmp{1}(ii, 1).Max = [];
    tmp{1}(ii, 1).DimensionsMode = 'Fixed';
    tmp{1}(ii, 1).SamplingMode = 'Sample based';    
    tmp{1}(ii, 1).SampleTime = -1;
    tmp{1}(ii, 1).DocUnits = '';
    tmp{1}(ii, 1).Description = ''; 

    outBus = Simulink.Bus;
    outBus.Description = '';
    outBus.DataScope = 'Auto';
    outBus.HeaderFile = '';
    outBus.Alignment = -1;
    outBus.Elements = tmp{1};    
    clear tmp;  

    %% Mem bus = outBus
    memBus = Simulink.Bus;
    memBus.Description = '';
    memBus.DataScope = 'Auto';
    memBus.HeaderFile = '';
    memBus.Alignment = -1;
    memBus.Elements = outBus.Elements;

    %% TP bus
    evalin('base','clear slBus1');
    TP = get_fmag_tunable_params;
    Simulink.Bus.createObject(TP);
    tpBus = evalin('base','slBus1');

    %% Final aggregation
    name = {};
    name{end+1} = 'SCDalgo_liuqenn_inBus';
    name{end+1} = 'SCDalgo_liuqenn_outBus';
    name{end+1} = 'SCDalgo_liuqenn_memBus';
    name{end+1} = 'SCDalgo_liuqenn_tpBus';

    bus = {};
    bus{end+1} = inBus;
    bus{end+1} = outBus;
    bus{end+1} = memBus;
    bus{end+1} = tpBus;

end