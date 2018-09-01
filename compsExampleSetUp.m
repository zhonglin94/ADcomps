function [ state, model, schedule ] = compsExampleSetUp
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%  Email: Yzhonglin_edu@163.com
% We define a 6x6x3 grid, spanning a 120*120*12 domain with homogenous
% permeability and porosity. The effect of gravity is ignored. For
% illustrative purpose, we added a simple well with a bhp control of 30
% bars and no well equations are needed and a methane injection with a gas
% rate of 1000m3/d
%
% SYNOPSIS:
%    [state, model, schedule] = compsExampleSetUp
%
% RETURNS:
%     state - The state of reservoir at a time step
%   model - Compositional model
% schedule - Include well and control information
%}
%% User Input
%     cartDims = [20,1,1];
%     physDims = [250, 100, 50].*0.6;
%     G = cartGrid(cartDims, physDims);
%     G = computeGeometry(G);
%     gravity reset off
%     
%     poro = 0.2;
%     perm = 20*milli*darcy;
%     rock = makeRock(G, perm, poro);
%% The SPE1 Example
  [G, rock] = setupSPE1();
    %% Input properties of each component
     %   CO2  C1      n-C4   n-C10
     zi =   [0.01 0.19 0.2 0.6]; %Mole Fraction    
     % Id represent the index of components in 'Components_db'
     id = [3 5 9 16];
     [mix, names] = getCompsProps(id);
     mix.x = zi;
    %% Set up wells and schedule
     T = 14000*day;
    % Producer
    W = [];    
    % Two wells are added to the model
     W = verticalWell(W, G, rock, 1, 1, (1:1), 'Name', 'P1', 'radius',  0.0762, ...
          'Type', 'bhp', 'Val', 140*barsa, 'Sign', -1);   % Producer
      W = verticalWell(W, G, rock, 10, 10, (3:3), 'Name', 'J1', 'radius',  0.0762, ...
          'Type', 'grat', 'Val', 300000*meter^3/day, 'Sign', 1, 'Comp_i', [0, 1, 0, 0, 0]); % Injector (CO2)
 %% run
     dT_target = 800*day; 
    dt = rampupTimesteps(T, dT_target, 10);
    % Set up a schedule
     schedule = struct();     
     schedule.step.val = dt; 
     schedule.step.control = (cumsum(dt) > T) + 1;
     schedule.control = struct('W', W);
     %%
     % Relative permeability
      props = getRelPermTable();     
     props.DENSITY = [0 1.0378e+03  0];
     props.PVTW = [2.768e+07,1.029,4.539e-10,3.100e-04,0];
     props.ROCK = [7.865e+02,1.038e+03,0.9698];
       props.PVTW = [2.768e+07,1.029,0 ,3.100e-04,0];
       props.ROCK = [7.865e+02,0 ,0.9698];
      fluid = initCompsADIFluid(props);  
     T = 60+273.15; % Reservoir temperature
     EoS = cPREoS; % Equation of state
     model = compositionalModel(G, rock, fluid, 'EOS', EoS, 'compsVarNames', names, 'T', T, 'Mixture', mix, 'Verbose', true); 

     % Water and hydrocarbon properties
     model.inputdata.props = props;
     model.inputdata.sinPhR = 1;
     model.inputdata.sinPhW = 2;
     %% Initial the state before simulation
     state = initResSol(G, 140*barsa, [0.4, 0, 0]);
     % Overall hydrocarbon composition
     state.zc = repmat(zi, model.G.cells.num,1); 
     state.ac = zeros(model.G.cells.num,1);
     % Implement flash for state and then get water, oil and gas properties
     state = flashState(model, state);    
end