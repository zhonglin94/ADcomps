% Required modules
mrstModule add ad-fi ad-props
% Setup 10x10x10 grid of 200x200x25 m model.
nx = 10;    ny = 10;    nz = 5; 
Dx = 200;   Dy = 200;   Dz = 25;
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G); 
% Assume homogeneous rock.
permX = 30*milli*darcy; poro  = 0.3;
rock.perm = repmat(permX, [G.cells.num, 1]);
rock.poro = repmat(poro , [G.cells.num, 1]);
% Set rock compressibility:
cr = 1e-6/barsa;
% The reference pore volume pv_r at the reference pressure p_r
pv_r = poreVolume(G, rock);
p_r  = 200*barsa;
% Finally, the pressure dependent function of pore-volumes
pv   = @(p) pv_r .* exp( cr .* (p - p_r) );
%% Fluid (oil) properties
mu   = 5*centi*poise; % Assume constant viscosity:
c    = 1e-3/barsa; % Assume constant oil compressibility
rho_r = 850*kilogram/meter^3; % Reference fluid density  at reference pressure p_r
rho   = @(p) rho_r .* exp( c .* (p - p_r) ); % The pressure dependent function of fluid density
rhoS = 750*kilogram/meter^3; % The fluid density at surface condition={1barsa, 298K}
%% Single horizontal well in J direction 
W=[]; W = verticalWell(W, G, rock, 1, 1, (1:5), 'Name', 'producer', 'radius',  0.0762, ...
         'Type', 'bhp', 'Val', 1*barsa, 'Sign', -1); % We consider a well with 5 connections in I direction. 
 %% Plotting
f = [figure(1), figure(2)];
set(0, 'CurrentFigure', f(1));
clf
plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
plotWell(G, W);
axis off;
set(f(1), 'Color', 'w'); 
camproj perspective;
view(3);
%% Initial conditions
% Ignore the gravity. Assume all grids share a same initial pressure
p_init = repmat(200*barsa, [G.cells.num, 1]);
%% Setting up components needed for the simulation
% Construct interior face connections (assume no-flow boundary condition)
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
% Two-point flux approximation with harmonic average to obtain the two-sided transmissibilities
hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);
% 'Gradient matrix' C as follows
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
% The discrete gradient and divergence operators for equations are now given by
grad = @(x)-C*x;
div  = @(x)C'*x;
% Average of cell-based quantities in neighboring cells
avg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
%% Pressure and well equations:
z = G.cells.centroids(:,3); % z-coordinate of grid cells
pressureEq = @(p, p0, dt) (1/dt) .* (pv(p).*rho(p) - pv(p0).*rho(p0)) ...
    - div( avg(rho(p) ./ mu) .* T .* grad(p) );
% Wellrates are given as Peaceman well-index times pressure drop
wc = W(1).cells; % perforation grid cells
WI = W(1).WI;    % well indices
dz = W(1).dZ;    % perforation depth relative to well reference depth
wellRates = @(p, bhp) WI .* (rho(p(wc)) ./ mu) .* (bhp - p(wc));
%% Define ADI variables
% The primary variables are grid-cell pressures, well bhp and surface rates
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
% For convenience, create indices to variables when stacked:
pIx = 1:G.cells.num; bhpIx = G.cells.num + 1; qSIx = G.cells.num + 2;
%% Set up simulation parameters
numSteps = 52; totTime  = 365*day; dt = totTime / numSteps;
% Tolerance and maximum number of iterations for the Newton solver.
tol      = 1e-5; maxits   = 10;

% Save output in array 'sol'
sol = repmat(struct('time', [], 'pressure', [], 'bhp', [], 'qS', []), ...
             [numSteps + 1, 1]);
sol(1).time     = 0;
sol(1).pressure = double(p_ad);
sol(1).bhp      = double(bhp_ad);
sol(1).qS       = double(qS_ad);

% Set up plot
set(0, 'CurrentFigure', f(2));
clf
set(f(2), 'Color', 'w')
subplot(2,1,1),
plotCellData(G, convertTo(p_init, barsa));
title('pressure [bar]', 'EdgeColor', 'w');
colorbar; view(3);
camproj perspective

subplot(2,1,2);
axis([0, convertTo(totTime,day), 0, 300]);
title('Surface volume rate [m^3/day]'); 
hold on
%% Main simulation
% Sloving the equations in a fully implicit manner with the help of
% automatic differentiation frameworks to compute the Jacobina matrix
t = 0;  step = 0; nDigits = floor(log10(maxits)) + 1;
while t < totTime,
    t = t + dt; step = step + 1;    
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
            step, convertTo(t - dt, day), convertTo(t, day));
    % Newton loop
    resNorm = 1e99;  p0  = double(p_ad); nit = 0;
    while (resNorm > tol) && (nit < maxits)
        % Create equations:
        eqs = cell([3, 1]);
        eqs{1} = pressureEq(p_ad, p0, dt);        
        eqs{1}(wc) = eqs{1}(wc) - wellRates(p_ad, bhp_ad); %  Constraint of well contributions
        eqs{2} = qS_ad - sum(wellRates(p_ad, bhp_ad))/rhoS; % Constraint of  wellrates
        eqs{3} = bhp_ad - 100*barsa; % Constraint of well bottom-hole pressure
        % Concatenate equations and solve:
        eq  = cat(eqs{:});
        J   = eq.jac{1};  % Jacobian
        res = eq.val;     % residual
        upd = -(J \ res); % Newton update
        % Update variables
        p_ad.val   = p_ad.val   + upd(pIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val  = qS_ad.val  + upd(qSIx);
 
        resNorm = norm(res);
        nit     = nit + 1;
        fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
    end
    
    if nit > maxits,
        error('Newton solves did not converge')
    else
        sol(step+1).time     = t;
        sol(step+1).pressure = double(p_ad);
        sol(step+1).bhp      = double(bhp_ad);
        sol(step+1).qS       = double(qS_ad);

        % Plot evolution
        set(0, 'CurrentFigure', f(2));
        subplot(2,1,1), cla, caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');

        subplot(2,1,2)
        plot(convertTo(sol(step+1).time, day), ...
             convertTo(-sol(step+1).qS , meter^3/day), '*');

        drawnow
    end
end
