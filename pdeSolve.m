function [ results, U ] = pdeSolve( model, tMin, tMax, plotStr, yLim )
% Solve advection diffusion problem
% u_t + au + bu_x = cu_xx + F
% du + eu_x = L, x = 0
% fu + gu_x = R, x = 1

% Impose initial condition
xAxis = model.grid;
U = model.initialCondition(xAxis);

% Solve PDE using ode45
RHS = model.RHS;
if nargin >= 4 && strcmp(plotStr,'plot')
    options = odeset('OutputFcn',@plotSol,'AbsTol',1e-10,'RelTol',1e-8);
    [~,U] = ode45(RHS,[tMin,tMax],U,options);
else
    options = odeset('AbsTol',1e-10,'RelTol',1e-8);
    [~,U] = ode45(RHS,[tMin,tMax/2,tMax],U,options);
end

results.grid = xAxis;
results.solution = U(end,:)';

% Exact solution and error
if ~isempty(model.solution)
    P = model.norm;
    solFcn = model.solution;
    solution = solFcn(xAxis,tMax);
    results.error = sqrt((U(end,:) - solution')*P*(U(end,:)' - solution));
else
    if nargout ~= 0
        warning('Analytic solution has not been provided')
    end
end

% Nested plot command
    function [ status ] = plotSol( t, y, flag )
        
        if strcmp(flag,'init')
            figure
            hold on
            set(gca,'NextPlot','replacechildren');
            set(gca,'YLim',yLim);
            set(gca,'XLim',[min(xAxis) max(xAxis)]);
            xlabel('$x$')
            ylabel('$u(x,t)$')
        end
        if ~strcmp(flag,'done')
            plot(xAxis,y,'k-');
            %
            %
            %axis([min(xAxis) max(xAxis) -0.5 2.5])
            drawnow;
        end
        status = 0;
        
    end

end