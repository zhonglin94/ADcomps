function [eqs, cq_mass, cq_s, mix_s, status, cstatus, cq_vol, cq_m] = computeWellContributionsSingleWellComps(wellmodel, wellSol, resmodel, q_s, pBH, packed)
[p, mob, ~ , dissolved] = unpackPerforationProperties(packed);

W = wellmodel.W;
assert(numel(wellSol) == 1);
assert(numel(W) == 1);
numPh = numel(q_s);
surfaceCondition = {1.*atm, (15.56+273.15).*Kelvin};

Tw = W.WI;
compi = W.compi;
wellStatus = W.status;
% Perforations/completions are closed if the well are closed or they are
% individually closed
perfStatus = W.cstatus.*wellStatus;
% Closed shut connection by setting WI = 0
Tw(~perfStatus) = 0;

% Well total volume rate at std conds:
qt_s = 0;
for ph = 1:numPh
    qt_s = qt_s + q_s{ph};
end
qt_s = qt_s*wellStatus;

% Get well signs, default should be that wells are not allowed to change sign
% (i.e., prod <-> inj)
if ~wellmodel.allowSignChange % injector <=> w.sign>0, prod <=> w.sign<0
    isInj = W.sign>0;
else
    isInj = double(qt_s)>0;   % sign determined from solution
end

%--------------------------------------------------------------------------
% Pressure drawdown (also used to determine direction of flow)
drawdown    = -(pBH+vertcat(wellSol.cdp)) + p;
connInjInx  = drawdown < 0; %current injecting connections

% A cross-flow connection is is defined as a connection which has opposite
% flow-direction to the well total flow
crossFlowConns = connInjInx ~= isInj;
% If crossflow is not alowed, close connections by setting WI=0
closedConns = ~wellSol.cstatus;
%closedConns    = false(size(crossFlowConns));
if ~wellmodel.allowCrossflow
    closedConns     = or(closedConns, crossFlowConns);
end
Tw(closedConns) = 0;
% Remove closedConns from connInjInx
connInjInx      = and(connInjInx, ~closedConns);

% ------------------ HANDLE FLOW INTO WELLBORE -------------------------
% producing connections phase volumerates:
cq_p = cell(1, numPh);
conEff = ~connInjInx.*Tw; % ok
for ph = 1:numPh
    cq_p{ph} = -conEff.*mob{ph}.*drawdown;
end
nc = resmodel.Mixture.numC;
xc = dissolved{1}; yc = dissolved{2}; rhom = dissolved{3};
cq_pm = cell(1, nc+1);
cqwm =  -conEff.*mob{1}.*drawdown.*rhom{1};
cq_pm{1} = cqwm;

[cqom, cqgm] = deal(cell(1,nc)); % Oil, gas molar rate
for i =1: nc
    cqom{i} = -conEff.*mob{2}.*drawdown.*rhom{2}.*xc{i};
    cqgm{i} = -conEff.*mob{3}.*drawdown.*rhom{3}.*yc{i};
    cq_pm{i+1} = cqom{i} +cqgm{i};
end

if ~all(connInjInx)
        [cq_ps, rhoS] = getWellRate(resmodel, cq_pm, connInjInx, surfaceCondition{:});
%         [cq_ps, rhoS]  = getCondtionedVolumeRate(resmodel, cq_pm, connInjInx, surfaceCondition);
else
    rhoS = cell(1,3);
    cq_ps = cell(1,3);
    for i = 1:numPh
        cq_ps{i} = 0;
        rhoS{i} = 0;
    end
end
%  test = cat(cq_ps{:});
% val = test.val
% jac = full(test.jac{:})
 
% producing connections phase volumerates at standard conditions:
% cq_ps = conn2surf(cq_p, b, dissolved, resmodel);
% Sum of phase rates from producing connections at std conds:
q_ps = cell(1, numPh);
for ph = 1:numPh
    q_ps{ph} = sum(cq_ps{ph});
end
%  double(cat(q_ps{:}))

isInj = double(qt_s)>0;
wbq = cell(1, numPh);
if isInj
    % Injection given by prescribed composition
    for ph = 1:numPh
        wbq{ph} = compi(ph).*qt_s;
    end
else
    % Determined by reservoir conditions
    for ph = 1:numPh
        wbq{ph} = q_s{ph}.*(q_s{ph}>0) - q_ps{ph};
    end
end

% compute wellbore total volumetric rates at std conds.
wbqt = 0;
for ph = 1:numPh
    wbqt = wbqt + wbq{ph};
end
% check for "dead wells":
deadWells = double(wbqt)==0;
if any(deadWells)
    for ph = 1:numPh
        wbq{ph} = wbq{ph}.*(~deadWells) + compi(ph).*deadWells;
        % Avoid division by zero
    end
    wbqt(deadWells) = 1;
end
% compute wellbore mixture at std conds
mix_s = cell(1, numPh);
for ph = 1:numPh
%     mix_s{ph} = wbq{ph}./wbqt;
        mix_s{ph} = 1;
end
% ------------------ HANDLE FLOW OUT FROM WELLBORE -----------------------
% total mobilities:
mt = mob{1};
for ph = 2:numPh
    mt = mt + mob{ph};
end
% injecting connections total volumerates
cqt_i = -(connInjInx.*Tw).*(mt.*drawdown);
% injecting connections total volumerates at standard condintions
if isInj && (wellSol.sign == 1)
        allcompi = W.compi;

        bhCondition = {pBH, resmodel.T};
        [compi, rhom] = getWellSat(resmodel, allcompi, bhCondition{:});
        
        
        cqt_im =(compi{1}.*rhom{1} + compi{2}.*rhom{2} +compi{3}.*rhom{3}).*cqt_i;
        compi = [double(compi{1}), double(compi{2}), double(compi{3})];
        cq_im = cell(1, nc+1);
        for  i = 1 : nc+1
           cq_im{i} = cqt_im.*allcompi(i); 
        end
        % injecting connections total volumerates at standard condintions
        cq_is = getWellRate(resmodel, cq_im, connInjInx, surfaceCondition{:});
%         [cq_is, ~] = getCondtionedVolumeRate(resmodel, cq_im, ~connInjInx, surfaceCondition);        
else
         cq_is = cell(1,3);
        for i = 1:numPh
            cq_is{i} = 0;
        end 
        cq_im = cell(1, nc+1);
        for i = 1 : nc+1
            cq_im{i} =0;
        end
end
% % connection phase volumerates at standard conditions (for output):
cq_s = cell(1,numPh);
for ph = 1:numPh
        cq_s{ph} = cq_ps{ph} + cq_is{ph};
end
cq_m = cell(1, nc+1);
for i = 1 : nc+1
    cq_m{i} = cq_pm{i} + cq_im{i};    
end

% Reservoir condition fluxes
cq_vol = cell(1, numPh);
for ph = 1:numPh
     cq_vol{ph} = connInjInx.*cqt_i.*compi(ph) + ~connInjInx.*cq_p{ph};
end
%---------------------- WELL EQUATIONS     -------------------------------
% Well equations
eqs = cell(1, numPh);
for ph = 1:numPh
    eqs{ph} = q_s{ph} - sum(cq_s{ph});
end
if ~all(wellStatus)
    % Overwrite equations with trivial equations for inactive wells
    subs = ~wellStatus;
    for ph = 1:numPh
        eqs{ph}(subs) = q_s{ph}(subs) - double(q_s{ph}(subs));
    end
end
% return mix_s(just values), connection and well status:
mix_s   = cell2mat( cellfun(@double, mix_s, 'UniformOutput', false));
cstatus = ~closedConns;
% For now, don't change status here
status = vertcat(wellSol.status);

cq_mass = cell(1, numPh);
for i = 1:numPh
       cq_mass{i} = cq_s{i}.*rhoS{i};
end

end
%--------------------------------------------------------------------------
%                     PRIVATE FUNCTION
function [compi, rhom] = getWellSat(model, allx, P, T)
    waterInj = true;    
    if ~all(allx(2:end) == 0)
        x= allx(2:end);
        x = x./sum(x);
        waterInj = false;
    end   
    
    f = model.fluid;
    pW =P; % Ignore capillary pressure
    rho1 = f.bW(pW).*f.rhoWS;                
    % Water, oil and gas weight
    mwW = 18.*milli; 
    rhom1 = rho1./mwW;   
    
   compi = {1, 0, 0};
   rhom = {rhom1, 0, 0};
   
 if ~waterInj
     [rhom, rho, beta] = flashWell(model, x, P, T);
     compi ={0, 1-beta, beta};
 end
end

%--------------------------------------------------------------------------

% function [compi, rhom, alpha] = getConditionedSaturation(resmodel, allcompi, condition)
%        allcompi = allcompi./sum(allcompi);
%         if sum(allcompi(2:end))>0
%             [~, rhom, alpha, ~]= flashWell(resmodel, allcompi(2:end), condition{:});
%         else
%             [rhoW, mwW] = getWaterDensity(resmodel, condition{1});
%             rhom ={rhoW/mwW, 1, 1}; alpha = 0.5; % Actually a water injection well
%         end
%         compW = allcompi(1);
%         compH = 1 - compW;
%         if isnumeric(rhom)
%             compi = [compW, (1-alpha).*compH, alpha.*compH]./rhom;
%             compi =compi ./ sum(compi);
%         else
%             compi = {compW, (1-alpha).*compH, alpha.*compH};
%             for i = 1 : 3
%                 compi{i} = compi{i}./ rhom{i};
%             end
%             compi = resmodel.normalize2Unit(compi);
%         end
% end
%---------------------------------------------------%
function [cqs, rho] = getWellRate(model, cqm, connInjInx, P, T)

    connNum = numel(connInjInx);
    cqs = cell(1, 3);
    rho = cell(1, 3);
    
    P = repmat(P, [connNum, 1]);
    
    f = model.fluid;
    pW =P; % Ignore capillary pressure
    rho1 = f.bW(pW).*f.rhoWS;                
    % Water, oil and gas weight
    mwW = 18.*milli; 
    rhom1 = rho1./mwW;   
    
    cqs{1}  = cqm{1}./rhom1;
    cqs{2}  =zeros(size(double(cqs{1})));
    cqs{3}  = zeros(size(double(cqs{1})));
    
    rho{1} = rho1;
    rho{2} = zeros(size(double(rho{1})));
    rho{3} = zeros(size(double(rho{1})));
    
    needFlash = false;
    temp = 0;
    for i = 2 : model.Mixture.numC+1
        temp = temp + cqm{i};
        if any(abs(double(cqm{i}))>0)
            needFlash = true;
        end
    end
    
    if needFlash
        x = cqm(2:end);
        [rhom, rho, beta] = flashWell(model, x, P, T); 

        mol2 = temp.*(1-beta);
        mol3 = temp.*beta;
        
        cqs{2} = mol2./rhom{2};
        cqs{3} = mol3./rhom{3};
    end
%     for i = 1 : 3
%         cqs{i}(connInjInx) = 0;
%         rho{i}(connInjInx) = 0;
%     end

end





