classdef compositionalModel < ReservoirModel
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% DESCRIPTION
% Full implicit compositional model
%}
 properties
     % Equation of state
     Eos    
     % Compositional variable names
     compsVarNames
     % Temperature of reservoir fluids
     T
     % Mixture of reservoir fluids
     Mixture
    % Maximum relative compositional variable increment
    dcMaxRel
     % Maximum relative compositional variable increment
    dcMaxAbs   
end   
methods
    function model = compositionalModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});
        % Equation of state
        model.Eos = [];
        model.compsVarNames =[];

        % Max increments of components
        model.dcMaxAbs = 0.2;
        model.dcMaxRel = 0.2;
        
        % Max increments of pressure
         model.dpMaxRel = 0.2;
%         model.dpMaxAbs = inf;
        % All phases are present 
        model.oil = true;
        model.gas = true;
        model.water = true;
        model.saturationVarNames = {'sw', 'so', 'sg'};       
        
        % Compositional -> use CNV style convergence 
        model.useCNVConvergence = true;         
        
        model = merge_options(model, varargin{:});
    end
    
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
%         [problem, state] = equationsCompositional(state0, state, model, dt, ...
%                         drivingForces, varargin{:}); 
        [problem, state] = equationsCompositional_WS(state0, state, model, dt, ...
                        drivingForces, varargin{:}); 
    end
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        compsVars = lower(model.compsVarNames);
        satVars = model.saturationVarNames;
        switch (lower(name))
            case compsVars
                    fn = 'zc';
                    index = model.getCompsVarsIndex(name);      
            case {'zc','xc','yc', 'ki'}
                fn = lower(name);
                index = 1:numel(compsVars);
            case {'alpha','rhoo','rhog', 'mwo','mwg','status', 'ac' }
                fn = lower(name);
                index = 1;         
            case{'rho','rhom', 'mw', 'mu'}
                fn = lower(name);
                index = 1:numel(satVars);
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
   end
   % --------------------------------------------------------------------%        
    function compsIndex = getCompsVarsIndex(model, name)
        compsVars = lower(model.compsVarNames);
        compsIndex = find(strcmpi(compsVars, name), true);
    end
    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
        ng = model.G.cells.num; % Number of grid
        nc = numel(model.compsVarNames); % Numer of components
            % zc must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'zc', [ng, nc], [1, 2]);
    end     
        % --------------------------------------------------------------------%    
        function state = flashState(model, state)
            sinPh = model.inputdata.sinPhR;
            T =  model.T; 
            P = state.pressure;            
            numG = model.G.cells.num; 
            
            mix = model.Mixture;
            xall = state.zc;            
            EoS = model.Eos;         
            fluid = model.fluid;
            % Flash for all grids

            [beta, xc, yc, ZL, ZG, st, findex] = flashALL(EoS, P, T, mix, xall,  numG, sinPh);

            state.xc = xc;
            state.yc = yc;            
            state.status = st;
            state.beta =beta;
            state.Z = [ZL, ZG];            
            
%             [rho, rhom, mw, mu] = updateFlash(model, P, T, state.Z, x, y);     
            [rhom, rho] = getDensity(mix, xc, yc, ZL, ZG, P, T, fluid);
            mu = getViscosity(mix, rhom, xc, yc, P, T, fluid);

            state.rho = rho;
            state.rhom =rhom;
            state.mu =mu;
            s1 = state.s(:,1);
             
             [s2, s3] = getOGSat(rhom, beta); % s2+s3 = 1
             
             s2 = s2.*(1-s1);
             s3 = s3.*(1-s1);
             s = [s1, s2, s3];            

             state.s = s;
        end
        % --------------------------------------------------------------------%
        function [rhom, rho, beta, st]= flashWell(model, x, P, T)   
            sinPh = model.inputdata.sinPhW;
            typeCell = true;
            if ~iscell(x)
                x = num2cell(x,1);
                 typeCell = false;
            end
              
            mix = model.Mixture;
            fluid = model.fluid;            
            
            numC = mix.numC;      
            numF = numel(double(x{1}));  % Number of connnections
            EoS = model.Eos;  
             
            x = normalize2unit(model, x, numC);
            xall = zeros(numF, numC);
            for i = 1 : numC
                xall(:, i) = double(x{i});                
            end
            
            P = repmat(P(1), [numF, 1]);           
            [beta, xc, yc, ZL, ZG, st, findex] = flashALL(EoS,  double(P), T, mix, xall, numF, sinPh);            
            [rhom, rho] = getDensity(mix, xc, yc, ZL, ZG, P, T, fluid);
            
            if ~iscell(rhom)
                rhom = num2cell(rhom, 1);
                rho = num2cell(rho, 1);
            end
        end        
        % --------------------------------------------------------------------%
        function x = normalize2unit(model, x, numC)
            temp = 0;
            for i  = 1 : numC
                temp = temp + x{i};
            end
            
            for i  = 1 : numC
                x{i} = x{i}./temp;
            end
        end
        % --------------------------------------------------------------------%  
        function [xcVars, ycVars] = splitPrimaryVariablesComps(model, vars)
            % Split a set of primary variables into three groups:
            % Well variables, saturation variables, gas composition variables and 
            % oil compOsition variables and the rest. This is useful because the
            % saturation variables oil and gas composition variables usually are updated
            % together, and the well variables are a special case.
            gvar = model.getCompsVarNames;
            xcVars = cellfun(@(x) strcat('x',x), gvar,'UniformOutput' , false);
            ycVars = cellfun(@(x) strcat('y',x), gvar,'UniformOutput' , false);

            isXc = cellfun(@(x) any(strcmpi(xcVars, x)), vars);
            isYc = cellfun(@(x) any(strcmpi(ycVars, x)), vars);
            xcVars = vars(isXc);
            ycVars = vars(isYc);
            for i = 1:numel(xcVars)
                vars = model.stripVars(vars,xcVars{i});
                vars = model.stripVars(vars,ycVars{i});
            end
        end   
        % --------------------------------------------------------------------%
        function vars = getCompsVarNames(model)
                    vars = model.compsVarNames;
        end
        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
             % Get the density and viscosity for oil and gas
            p =state.pressure;
            numG = model.G.cells.num;
            xc =state.xc;
            yc =state.yc;
            state.ki = xc./yc;
            
            EoS = model.Eos;     
            T =  model.T; 
            mix = model.Mixture;
            fluid = model.fluid;
            
             [fL, fugL, ZL] = fugF(EoS, T, p, mix, xc , 'liq'); 
             [fG, fugG, ZG] = fugF(EoS, T, p, mix, yc , 'gas'); 
             [rhom, rho] = getDensity(mix, xc, yc, ZL, ZG, p, T, fluid);
             mu = getViscosity(mix, rhom, xc, yc, p, T, fluid);

            state.rho = rho;
            state.rhom = rhom;
            state.mu = mu;             
            
            s2 = state.s(:, 2);
            s3 = state.s(:, 3);
            rhom2 = rhom(:, 2);
            rhom3 = rhom(:, 3);
            mol2 = rhom2.*s2;
            mol3 = rhom3.*s3;
            molt = mol2+mol3;
            
            state.beta = mol3./molt;
            
            for i = 1 : numG
                state.zc(i, :) = (mol2(i).*xc(i, :) + mol3(i).*yc(i, :))./molt(i);
            end
            state = getTotalMolar(model,state);
             stbb = state.status; 
            if ~all(state.status == stbb)
                stop = true;
            end
            report =[];
        end
          % --------------------------------------------------------------------%
         function changed = phaseChangedInWellBlocks(model, state, drivingForces)   
             nw = numel(drivingForces.W);  % Number of well
             for i = 2 : nw
                 nwc = numel(drivingForces.W(i).cells); % Number of connnection in a well
                 wName = drivingForces.W(i).name;
                 for j = 1 : nwc
                     wc = drivingForces.W(i).cells(j);
                     zc = state.zc(wc, :); p = state.pressure(wc);
                     try
                         [ ~, ~, stw, ~]  = flasher(model, zc, p);
                     catch
                         error('Error stability check for block %d in well %s', wc, wName)
                     end
                     if ~(state.status(wc) == stw)
                         changed = true;
                         return
                     end
                 end
             end
             changed = false;
         end
        % --------------------------------------------------------------------%
        function  output = assignValue(model, input0, input, identifier)
            input0(identifier, :) = input(identifier, :);
            output = input0;      
        end
        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
         [xcVars, ycVars] = model.splitPrimaryVariablesComps(problem.primaryVariables);
         % Update mole fraction in one go
        state  = model.updateXc(state, dx, problem, xcVars);   
        state  = model.updateYc(state, dx, problem, ycVars);       
%         state = model.updateSat(state, dx, problem, satVars);
        
        vars = problem.primaryVariables;
        removed = false(size(vars));
        [vars, ix] = model.stripVars(vars, {xcVars{:}, ycVars{:}});
        removed(~removed) = removed(~removed) | ix;
        
        problem.primaryVariables = vars;
        dx(removed) = [];
         % Parent class handles almost everything for us
        [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
        end
        % --------------------------------------------------------------------%        
        function state = updateSat(model, state, dx, problem, satVars)
                ds = zeros(model.G.cells.num, numel(satVars)+1);
                for i =2:numel(satVars)+1
                    ds(:,i) = model.getIncrement(dx, problem, satVars{i-1}); % Water, oil and gas
                end
                ds(:,1) = -sum(ds(:,2:end),2);
                state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);
        end
    % --------------------------------------------------------------------%        
    function state = updateXc(model, state, dx, problem, xcVars)
        if nargin < 5
            % Get the xc names directly from the problem
            [xcVars] = ...
                splitPrimaryVariablesComps(model, problem.primaryVariables);
        end        
        if isempty(xcVars)
            % No saturations passed, nothing to do here.
            return
        end        
        dxc = zeros(model.G.cells.num, numel(xcVars)+1);
        for i = 1:numel(xcVars)
             v = model.getIncrement(dx, problem, xcVars{i});
             dxc(:,i) =v;
        end
         dxc(:,end) = -sum(dxc,2);
         state = model.updateStateFromIncrement(state, dxc, problem, 'xc', model.dcMaxRel, model.dcMaxAbs);
        % We update all mole fractions simultanously    
    end
    % --------------------------------------------------------------------%   
    function state = updateYc(model, state, dx, problem, ycVars)
        if nargin < 5
            % Get the yc names directly from the problem
            [~, ycVars] = ...
                splitPrimaryVariablesComps(model, problem.primaryVariables);
        end        
        if isempty(ycVars)
            % No saturations passed, nothing to do here.
            return
        end        
        dyc = zeros(model.G.cells.num, numel(ycVars)+1);
        for i = 1:numel(ycVars)
             v = model.getIncrement(dx, problem, ycVars{i});
             dyc(:,i) =v;
        end
         dyc(:,end) = -sum(dyc,2);
         state = model.updateStateFromIncrement(state, dyc, problem, 'yc', model.dcMaxRel, model.dcMaxAbs);
        % We update all mole fractions simultanously    
    end
    % --------------------------------------------------------------------%   
    function state = setProp(model, state, name, value)
        [fn, index] = model.getVariableField(name);
        if strcmpi(fn,'zc')
            state.(fn)(:, index) = num2cell(value,1);
        else
            state = setProp@PhysicalModel(model, state, name, value);
        end
    end
    % --------------------------------------------------------------------%
    function [flowEqsNames, fugEqsNames] = getEqNames(model)
        compsVars = model.getCompsVarNames;
        flowEqsNames = cellfun(@(compsVars) strcat('flow',compsVars), compsVars,'UniformOutput' , false);
        flowEqsNames = {'flowWater', flowEqsNames{:}};
        fugEqsNames =cellfun(@(compsVars) strcat('fug',compsVars), compsVars,'UniformOutput' , false);
    end
    % --------------------------------------------------------------------%

    function [eqs, names, types, wellSol, src] = insertWellEquations(model, eqs, names, srcEqNames,...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     dissolved, components, ...
                                                     dt, opt)
        % Add in the effect of wells to a system of equations, by adding
        % corresponding source terms and augmenting the systme with
        % additional equations for the wells.
        %
        % INPUT
        %
        % eqs    - Cell array of equations that are to be updated.
        %
        % names  - The names of the equations to be updated. If
        %          phase-pseudocomponents are to be used, the names must
        %          correspond to some combination of "water", "oil", "gas"
        %          if no special component treatment is to be introduced.
        %
        % types  - Cell array with the types of "eqs". Note that these
        %          types must be 'cell' where source terms is to be added.
        %
        % src    - Struct containing all the different source terms that
        %          were computed and added to the equations.
        %
        % Remaining input arguments correspond to a variety of reservoir
        % properties with self-explanatory names.
        if model.FacilityModel.getNumberOfWells() == 0
            return
        end
        fm = model.FacilityModel;
        nc = numel(model.getCompsVarNames);
        [src, wellsys, wellSol] = ...
            fm.getWellContributions(wellSol0, wellSol, wellVars, ...
                                    wellMap, p, mob, rho, dissolved, components, ...
                                    dt, opt.iteration);
                               
        wc = src.sourceCells;
        % Treat phase pseudocomponent source terms from wells
        for i = 1:nc+1
            sub = strcmpi(names, srcEqNames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                eqs{sub}(wc) = eqs{sub}(wc) - src.phaseMolar{i};
            end
        end
        % Treat component source terms from wells
        cnames = model.getComponentNames();
        for i = 1:numel(cnames)
            sub = strcmpi(names, cnames{i});
            if any(sub)
                assert(strcmpi(types{sub}, 'cell'), 'Unable to add source terms to equation that is not per cell.');
                eqs{sub}(wc) = eqs{sub}(wc) - src.components{i};
            end
        end
        offset = numel(wellsys.wellEquations);
        eqs(end+1:end+offset) = wellsys.wellEquations;
        names(end+1:end+offset) = wellsys.names;
        types(end+1:end+offset) = wellsys.types;
        eqs{end+1} = wellsys.controlEquation;
        names{end+1} = 'closureWells';
        types{end+1} = 'well';
    end
        % --------------------------------------------------------------------%
    function model = validateModel(model, varargin)
        if isempty(model.FacilityModel)
            model.FacilityModel = FacilityModelComps(model); %#ok
        end
        if nargin > 1
            W = varargin{1}.W;
            model.FacilityModel = model.FacilityModel.setupWells(W);
        end
        model = validateModel@PhysicalModel(model, varargin{:});
        return
    end

    % --------------------------------------------------------------------%
end
end

       