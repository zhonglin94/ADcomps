function [problem, state] = equationsCompositional_WS(state0, state, model, dt, drivingForces, varargin)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
% Email: Yzhonglin_edu@163.com
% DESCRIPTION
% Generate linearized problem for the compositional equations which 
% include the conservation equations, fugacity equations and well equations      
 %}

opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'iteration',   -1);
opt = merge_options(opt, varargin{:});
% Shorter names for some commonly used parts of the model and forces.
fluid = model.fluid;
o = model.operators;
nc = numel(model.compsVarNames);
numG = model.G.cells.num;
% 
gvar = model.compsVarNames;
xcVars = cellfun(@(gvar) strcat('x',gvar), gvar,'UniformOutput' , false);
ycVars = cellfun(@(gvar) strcat('y',gvar), gvar,'UniformOutput' , false);
% Properties at current timestep

[p, sO, sG, xc, yc, zc, ac, st, wellSol] = model.getProps(state, ...
    'pressure', 'oil', 'gas', 'xc', 'yc', 'zc', 'ac', 'status', 'wellSol');
% Properties at previous timestep
[p0, sO0, sG0, xc0, yc0, ac0, rhom0, wellSol0] = model.getProps(state0, ...
    'pressure', 'oil', 'gas', 'xc', 'yc', 'ac','rhom', 'wellSol');

xc0 = num2cell(xc0,1);
yc0 = num2cell(yc0,1);
s0 = {1-sO0-sG0, sO0,sG0};
rhom0 = num2cell(rhom0,1);

xc = num2cell(xc,1);
yc = num2cell(yc,1);
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize primary variabes, include fluid compositions
 [p, xc{1:end-1},sO, yc{1:end-1}, sG, ac, wellVars{:}] = initVariablesADI(p,  xc{1:end-1}, sO, yc{1:end-1}, sG, ac, wellVars{:}); 
tempx =0; tempy =0;
for i = 1:nc-1
    tempx =tempx+xc{:,i};
    tempy =tempy+yc{:,i};
end
    xc{:,end} = 1-tempx;
    yc{:,end} = 1-tempy;
    s = {1-sO-sG, sO, sG};
% Get CO2 fugacity in aqueous phase
 aw = 1-ac;
 aw0 = 1-ac0;
 fugC = CO2fugacity(model, ac, p);
 %
 T = model.T;
%  fugC{1} = fugC{1}.*kSalt(T);

 % Get the gas and oil properties and derivatives for ADI variables
 %%
 EoS = model.Eos;     
 T =  model.T; 
 mix = model.Mixture;
 
 [fL, fugL, ZL] = fugF(EoS, T, p, mix, xc , 'liq'); 
 [fG, fugG, ZG] = fugF(EoS, T, p, mix, yc , 'gas'); 
 
[rhom, rho] = getDensity(mix, xc, yc, ZL, ZG, p, T, fluid);
mu = getViscosity(mix, rhom, xc, yc, p, T, fluid);
%% Phase disappearance
tol = 1e-4;
disIndex2 = s{2}<tol;
disIndex3 = s{3}<tol;

st(disIndex2) = 2;
st(disIndex3) = 1;

% Evaluate relative permeability
[krW, krO, krG] = model.evaluateRelPerm(s);
T = o.T;
%Water
mobW   = krW./mu{1};
pW =p;
dpW    = o.Grad(pW);
% water upstream-index
upcw  = (double(dpW)<=0);
mobWcomps = {mobW.*aw, mobW.*ac};
vW{1} = -o.faceUpstr(upcw, mobWcomps{1}).*T.*dpW;
vW{2} = -o.faceUpstr(upcw, mobWcomps{2}).*T.*dpW;
% Oil
dpO    = o.Grad(p);
mobO   = krO./mu{2};
% Oil upstream-index
upco = (double(dpO)<=0);
vO = cell(1, nc);
for i = 1:nc
    mobOcomps = mobO.*xc{i};
    vO{i}    = - o.faceUpstr(upco, mobOcomps).*T.*dpO;
end   
% Gas
pG = p;
mobG   = krG./mu{3};
dpG    = o.Grad(pG);
    % Gas upstream-index
upcg    = (double(dpG)<=0);
vG = cell(1,nc);
for i =1:nc
   mobGcomps = mobG.*yc{i};
   vG{i} = - o.faceUpstr(upcg, mobGcomps).*T.*dpG;
end
mob = {mobW, mobO, mobG};
%% Equations 
% The index 1 stands for water
% water = (o.pv/dt).*( rhom{1}.*s{1} - rhom0{1}.*s0{1}) + o.Div(o.faceUpstr(upcw, rhom{1}).*vW);
water = (o.pv/dt).*( rhom{1}.*s{1}.*aw - rhom0{1}.*s0{1}.*aw0) + o.Div(o.faceUpstr(upcw, rhom{1}).*vW{1});
aCO2 = (o.pv/dt).*( rhom{1}.*s{1}.*ac - rhom0{1}.*s0{1}.*ac0) + o.Div(o.faceUpstr(upcw, rhom{1}).*vW{2});

[hydrocarbon, fugacity] = deal(cell(nc,1));
for i =1:nc
% Molar balance equations
% The index 2 and 3 stand for oil and gas respectively
hydrocarbon{i} = (st == 3 ).*((o.pv/dt).*( rhom{2}.*s{2}.*xc{i} - rhom0{2}.*s0{2}.*xc0{i}) + o.Div(o.faceUpstr(upcw, rhom{2}).*vO{i})+...
                          (o.pv/dt).*( rhom{3}.*s{3}.*yc{i} - rhom0{3}.*s0{3}.*yc0{i} ) + o.Div(o.faceUpstr(upcw, rhom{3}).*vG{i})) + ...
                          (st == 2) .* ((o.pv/dt).*( rhom{3}.*s{3}.*yc{i} - rhom0{3}.*s0{3}.*yc0{i} ) + o.Div(o.faceUpstr(upcw, rhom{3}).*vG{i})) +...
                          (st == 1) .* ((o.pv/dt).*( rhom{2}.*s{2}.*xc{i} - rhom0{2}.*s0{2}.*xc0{i}) + o.Div(o.faceUpstr(upcw, rhom{2}).*vO{i}));
% Fugacity equations
% fugacity{i} = (st == 3 ).*(fugL{i}-fugG{i}) + ...
%                     (st == 2).*(xc{i} - zc(:, i)) + ... % A trivial equation
%                     (st == 1).*(yc{i} - zc(:, i));  


fugacity{i} = (st == 3 ).*(fugL{i}-fugG{i}) + ...
                    (st == 2).*(xc{i} - xc0{i}) + ... % A trivial equation
                    (st == 1).*(yc{i} - yc0{i});   

end
hydrocarbon{1} = hydrocarbon{1}  + aCO2;
fugacity{nc} = (st == 2) .* (s{2} -0) + (st == 1) .* (s{3} -0) +...
                      (st == 3 ).*(fugL{nc}-fugG{nc});          
fugacity{nc+1} =   fugL{1} - fugC{1};
%% Put the set of equations into cell arrays along with their names/types.
eqs = {water, hydrocarbon{:}, fugacity{:}};
names = {'pressure', xcVars{1:end-1}, 'so', ycVars{1:end-1}, 'sg','ac'};
srcEqNames = {'pressure', xcVars{1:end-1}, 'so'};
primaryVars = {names{:}, wellVarNames{:}};
ncell = 2.*nc+1 + 1;
types = cell(1, ncell);
% 2nc+1 variables
for i =1:ncell
types{i} = 'cell'; 
end
dissolved = {xc, yc, rhom}; 
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, srcEqNames, types,... 
                                                      wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
