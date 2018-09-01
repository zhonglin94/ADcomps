classdef FacilityModelComps < FacilityModel
    
    methods
        function model = FacilityModelComps(reservoirModel, varargin)
            model = model@FacilityModel(reservoirModel, varargin{:});
            
        end
        % -------------------------------------------------------%
       function [sources, wellSystem, wellSol] = getWellContributions(model, wellSol0, wellSol, wellvars, wellMap, p, mob, rho, dissolved, comp, dt, iteration)
            % Get the source terms due to the wells, control and well
            % equations and updated well sol. Main gateway for adding wells
            % to a set of equations.
            if isnan(iteration) || iteration < 0
                warning(['Iteration number is not passed on to well model,', ...
                         'this may indicate wellbore pressure-drop will never be updated']);
            end
            actWellIx = model.getIndicesOfActiveWells();
            nw = numel(actWellIx);

            % Stores base well equations describing reservoir coupling
            allBaseEqs = cell(nw, 1);
            % Control equations ensure that we enforce constraints
            allCtrl = cell(nw, 1);
            % Volumetric phase source terms
            allVol = cell(nw, 1);
            % Molar phase source terms
            allMol = cell(nw, 1);
            % Mass phase source terms
            allMass = cell(nw, 1);
            % Composition source terms
            allComp = cell(nw, 1);
            

            % Get the additional equations not implemented in the minimal
            % "SimpleWell" class.
            enames = model.addedEquationNames;
            etypes = model.addedEquationTypes;
            cnames = model.ReservoirModel.getComponentNames();
            ncomp = numel(cnames);
            assert(ncomp == numel(comp), ...
                'Number of input components must match length of getComponentNames!');

            n_extra = numel(enames);
            assert(numel(etypes) == n_extra);

            allExtraEqs = cell(nw, n_extra);
            resModel = model.ReservoirModel;

            addedVars = model.addedPrimaryVarNames;
            varmaps = cell(1, numel(addedVars));
            for varNo = 1:numel(addedVars)
                varmaps{varNo} = model.getWellVariableMap(addedVars{varNo});
            end
            
            isBH = wellMap.isBHP;
            isQ = wellMap.isRate;
            emap = wellMap.extraMap;
            
            bhp = wellvars{isBH};
            qWell = wellvars(isQ);
            wellvars = wellvars(~(isBH | isQ));

            [basenames, basetypes] = model.WellModels{1}.getWellEquationNames(resModel);
            for i = 1:nw
                wellNo = actWellIx(i);
                wm = model.WellModels{wellNo};
                ws = wellSol(wellNo);
                ws0 = wellSol0(wellNo);

                W = wm.W;
                packed = packPerforationProperties(W, p, mob, rho, dissolved, comp, wellvars, addedVars, varmaps, emap, i);
                qw = cellfun(@(x) x(i), qWell, 'uniformoutput', false);
                bh = bhp(i);
                % Update pressure
                % ws = wm.updateConnectionPressureDrop(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                % Update limits
                [qw, bh, ws, ok] = wm.updateLimits(ws0, ws, resModel, qw, bh, packed, dt, iteration);
                if ~ok
                    bhp(i) = bh;
                    for phNo = 1:numel(qw)
                        qWell{phNo}(i) = qw{phNo};
                    end
                end
               % Set up well equations and phase source terms
               [allBaseEqs{i}, allCtrl{i}, extraEqs, extraNames, allMass{i}, allVol{i}, allMol{i}, ws] =...
                   wm.computeWellEquationsComps(ws0, ws, resModel, qw, bh, packed, dt, iteration);

               % Get component source terms and corresponding equations (if
               % any components are present)
               [compEqs, allComp{i}, compNames, ws] =...
                   wm.computeComponentContributions(ws0, ws, resModel, qw, bh, packed, allMass{i}, allVol{i}, dt, iteration);

               extraEqs = {extraEqs{:}, compEqs{:}};
               extraNames = {extraNames{:}, compNames{:}};

               for eqNo = 1:numel(extraEqs)
                   % Map into global list of equations
                   ix = strcmpi(enames, extraNames{eqNo});
                   allExtraEqs{i, ix} = extraEqs{eqNo};
               end
               wellSol(wellNo) = ws;
            end
            % We have assembled all equations for each well. Combine the
            % equations from the different wells into one (array) of each
            % type.
            nPh = nnz(resModel.getActivePhases);
            [srcMass, srcVol, eqs] = deal(cell(1, nPh));
            for phNo = 1:nPh
                srcMass{phNo} = combineCellData(allMass, phNo);
                srcVol{phNo} = combineCellData(allVol, phNo);
                eqs{phNo} = combineCellData(allBaseEqs, phNo);
            end
            
            nc = numel(resModel.getCompsVarNames);
            srcMol = cell(1, nc+1);
            for ncNo = 1:nc+1 % include water
                srcMol{ncNo} =  combineCellData(allMol, ncNo);
            end
            % Components are ordered canonically by reservoir model
            srcComp = cell(1, ncomp);
            for cNo = 1:ncomp
                srcComp{cNo} = combineCellData(allComp, cNo);
            end
            % If we have extra equations, add them in
            extraEqs = cell(1, n_extra);
            for i = 1:n_extra
                ok = ~cellfun(@isempty, allExtraEqs(:, i));
                extraEqs{i} = vertcat(allExtraEqs{ok, i});
            end
            % Equations are the base, common variables as well as any extra
            % equations added due to complex wells.
            names = horzcat(basenames, enames);
            types = horzcat(basetypes, etypes);

            eqs = {eqs{:}, extraEqs{:}};
            ctrleq = vertcat(allCtrl{:});

            wc = model.getActiveWellCells();
            [wc, srcMass, srcVol] = model.handleRepeatedPerforatedcells(wc, srcMass, srcVol);
            wellSystem = struct('wellEquations', {eqs}, ...
                                'names',  {names}, ...
                                'types', {types}, ...
                                'controlEquation', ctrleq);
            sources = struct('phaseMass',   {srcMass}, ...
                             'phaseVolume', {srcVol}, ...
                             'phaseMolar', {srcMol},...
                             'components',  {srcComp},...
                             'sourceCells', wc);
            if model.ReservoirModel.extraWellSolOutput
                wellSol = model.setWellSolStatistics(wellSol, sources);
            end
       end
        %----------------------------------------------%
        function model = setupWells(model, W, wellmodels)
            % Set up well models for changed controls or the first
            % simulation step.
            %
            % INPUT:
            % 
            % W       - Well struct (obtained from e.g. addWell or
            %           processWells)
            %
            % wellmodels (OPTIONAL ARGUMENT)
            %          - Cell array of equal length to W, containing class
            %          instances for each well (e.g. SimpleWell,
            %          MultisegmentWell, or classes derived from these). 
            %          If not provided, well models be constructed from the
            %          input. 
            nw = numel(W);
            if model.getNumberOfWells == 0
                % First time setup
                [pvars, fromRes, eqnames, eqtypes] = deal(cell(nw, 1));
                model.WellModels = cell(nw, 1);
                for i = 1:nw
                    % Set up models. SimpleWell for the time being
                    if nargin < 3
                        if isfield(W(i), 'isMS') && W(i).isMS
                            wm = MultisegmentWellComps(W(i));
                        else
                            wm = SimpleWellComps(W(i));
                        end
                    else
                        wm = wellmodels{i};
                    end
                    wm.dsMaxAbs = model.ReservoirModel.dsMaxAbs;
                    wm.dpMaxAbs = model.ReservoirModel.dpMaxAbs;
                    wm.dpMaxRel = model.ReservoirModel.dpMaxRel;

                    if isfield(W(i), 'vfp_index')
                        vfp_ix = W(i).vfp_index;
                        if vfp_ix > 0
                            if wm.isInjector()
                                vfp = model.VFPTablesInjector{vfp_ix};
                            else
                                vfp = model.VFPTablesProducer{vfp_ix};
                            end
                            wm.VFPTable = vfp;
                        end
                    end
                    % Get the added primary variables for this well, plus
                    % the equations and equation types it adds
                    [pvars{i}, fromRes{i}] = wm.getExtraPrimaryVariableNames(model.ReservoirModel);
                    [eqnames{i}, eqtypes{i}] = wm.getExtraEquationNames(model.ReservoirModel);
                    model.WellModels{i} = wm;
                end
                % Combine the different equations and types added by the
                % different wells into a canonical ordering.
                [model.addedPrimaryVarNames, keepix] = uniqueStable([pvars{:}]);
                model.addedPrimaryVarNamesIsFromResModel = [fromRes{:}];
                model.addedPrimaryVarNamesIsFromResModel = model.addedPrimaryVarNamesIsFromResModel(keepix);
                
                [model.addedEquationNames, keepix] = uniqueStable([eqnames{:}]);

                etypes = [eqtypes{:}];
                model.addedEquationTypes = etypes(keepix);
            else
                assert(model.getNumberOfWells == nw, ...
                    'Number of wells in facility model has changed during simulation')
                for i = 1:nw
                    % Update with new wells. Typically just a quick
                    % overwrite of existing wells
                    model.WellModels{i} = model.WellModels{i}.updateWell(W(i));
                end
            end
        end
        %----------------------------------------------% 
    end
end

function d = combineCellData(data, ix)
    d = cellfun(@(x) x{ix}, data, 'UniformOutput', false);
    isAD = ~cellfun(@isnumeric, d);
    if any(~isAD) && any(isAD)
        s = d(isAD);
        s = s{1};
        for i = 1:numel(isAD)
            if ~isAD(i)
                d{i} = double2ADI(d{i}, s);
            end
        end
    end
    d = vertcat(d{:});
end
