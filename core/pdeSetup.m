classdef pdeSetup < handle
    %PDESETUP Input class for solving PDEs
    %   Input for advection diffusion problem
    %   u_t + au + bu_x = cu_xx + F
    %   du + eu_x = L, x = L
    %   fu + gu_x = R, x = R
    
    properties
        dimension
        geometry
        discretisation
        coefficient
        forcing
        boundaryCondition
        initialCondition
        interfaceCondition
        solution
    end
    
    properties (Hidden)
        RHS
        grid
        norm
    end
    
    methods
        % Constructor function
        function [ model ] = pdeSetup( d )
            if ~(d==1 || d==2)
                error('Only 1 or 2 dimensions supported')
            else
                model.dimension = d;
            end
        end
        
        % Geometric functions
        function initGeometry( model, d )
            if nargin == 1
                model.geometry = pdeGeometry(model.dimension);
            else
                model.geometry = pdeGeometry(d);
                model.dimension = d;
            end
        end
        
        function setBoundary( model, vertices, faces )
            addBoundary(model.geometry, vertices, faces);
        end
        
        function setInterface( model, vertices, faces )
            if model.dimension == 1
                addInterface(model.geometry, vertices);
            elseif model.dimension == 2
                addInterface(model.geometry, vertices, faces);
            else
                error('Dimension error: Only 1 or 2 dimensions supported')
            end
        end
        
        function setRegion( model, nodes )
            defineRegion(model.geometry, nodes);
        end
        
        function plotModel( model )
            plotGeometry( model.geometry )
        end
        
        % Discretisation functions
        function initDiscretisation( model, d )
            if isempty(model.geometry.region)
                if nargin == 1
                    model.discretisation = pdeDiscretisation(model.dimension);
                else
                    model.discretisation = pdeDiscretisation(d);
                    model.dimension = d;
                end
            else
                nRegions = length(fieldnames(model.geometry.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    if nargin == 1
                        model.discretisation.(fieldName) = ...
                            pdeDiscretisation(model.dimension);
                    else
                        model.discretisation.(fieldName) = ...
                            pdeDiscretisation(d);
                        model.dimension = d;
                    end
                end
            end
        end
        
        function setResolution( model, Nx )
            if model.dimension ~= 1
                error('Dimension error: Only 1 dimension supported')
            end
            if isa(model.discretisation,'pdeDiscretisation')
                setRes(model.discretisation,Nx);
            else
                nRegions = length(fieldnames(model.geometry.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    setRes(model.discretisation.(fieldName),Nx);
                end
            end
        end
        
        function setRegRes( model, region, Nx )
            setRes(model.discretisation.(region),Nx);
        end
        
        function setOperator( model, op, file, dir )
            if model.dimension ~= 1
                error('Dimension error: Only 1 dimension supported')
            end
            if nargin == 3
                dir = cd;
            end
            if isa(model.discretisation,'pdeDiscretisation')
                setOp(model.discretisation,op,file,dir);
            else
                nRegions = length(fieldnames(model.geometry.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    setOp(model.discretisation.(fieldName),op,file,dir);
                end
            end
        end
        
        function setRegOp( model, region, op, file, dir )
            if nargin == 4
                dir = cd;
            end
            setOp(model.discretisation.(region),op,file,dir);
        end
        
        function getStep( model )
            if model.dimension ~= 1
                error('Dimension error: Only 1 dimension supported')
            end
            if isa(model.discretisation,'pdeDiscretisation')
                calcStep(model.discretisation,...
                    model.geometry.boundary.vertices);
            else
                nRegions = length(fieldnames(model.geometry.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    calcStep(model.discretisation.(fieldName),...
                        model.geometry.region.(fieldName));
                end
            end
        end
        
        % PDE coefficients
        function setCoefficients( model, a, b, c )
            if ~(numel(a)==1 && numel(b)==1 && numel(c)==1)
                error('Input error: Only scalar inputs supported')
            end
            if model.dimension > 1
                error('Dimension error: Only 1 dimension supported')
            else
                model.coefficient = [a, b, c];
            end
        end
        
        % Forcing function
        function setForcing( model, F )
            if ~isa(F,'function_handle')
                error('Input error: Expected function handle input')
            end
            model.forcing = F;
        end
        
        % Boundary conditions
        function setBC( model, data, varargin )
            % Check the number of input parameters
            narginchk(4,6);
            
            % Count possible boundary nodes
            nNodes = size(model.geometry.boundary.vertices,1);
            
            % Validate input
            p = inputParser;
            expectedTypes = {'Dirichlet','Neumann','Robin'};
            expectedNodes = {'B1'};
            for k=1:nNodes
                expectedNodes{k} = ['B',num2str(k)];
            end
            
            addRequired(p,'model', @(x) isa(x,'pdeSetup'));
            addRequired(p,'data', @(x) isa(x,'function_handle'));
            addRequired(p,'node',...
                @(x) any(validatestring(x,expectedNodes)));
            addRequired(p,'penalty',@isnumeric);
            addRequired(p,'BCtype',...
                @(x) any(validatestring(x,expectedTypes)));
            addOptional(p,'coefficients',[0 0],...
                @(x) validateattributes(x,{'numeric'},{'numel',2,'nonzero'}));
            
            parse(p, model, data, varargin{:});
            
            if strcmp(p.Results.BCtype,'Robin') && ~any(p.Results.coefficients)
                error(['Input error: ',...
                    'Robin condition requires two non-zero coefficients'])
            end
            
            BCstruct = struct;
            BCstruct.data = p.Results.data;
            BCstruct.node = p.Results.node;
            BCstruct.penalty = p.Results.penalty;
            BCstruct.BCtype = p.Results.BCtype;
            if strcmp(BCstruct.BCtype,'Robin')
                BCstruct.coefficients = p.Results.coefficients;
            end
            
            if isempty(model.boundaryCondition)
                model.boundaryCondition.BC1 = BCstruct;
            else
                nBC = length(fieldnames(model.boundaryCondition));
                fieldName = ['BC',num2str(nBC+1)];
                model.boundaryCondition.(fieldName) = BCstruct;
            end
        end
        
        % Initial condition
        function setIC( model, fun )
            if ~isa(fun,'function_handle')
                error('Input error: Expected function handle input')
            end
            model.initialCondition = fun;
        end
        
        % Interface condition
        function setIF( model, node, sigma )
            if ~ischar(node) || node(1)~='I' || length(node)~=2
                error('Input error: Expected input of form In')
            end
            if numel(sigma) ~= 2
                error('Input error: Two penalties expected')
            end
            IFstruct = struct;
            IFstruct.node = node;
            IFstruct.penalty = sigma;
            if isempty(model.interfaceCondition)
                model.interfaceCondition.IF1 = IFstruct;
            else
                nIF = length(fieldnames(model.interfaceCondition));
                fieldName = ['IF',num2str(nIF+1)];
                model.interfaceCondition.(fieldName) = IFstruct;
            end
        end
        
        % Analytic solution
        function setSolution( model, sol )
            if ~isa(sol,'function_handle')
                error('Input error: Expected function handle input')
            end
            model.solution = sol;
        end
        
        % Prepare PDE for solution
        function pdePrepare( model )
            if model.dimension ~= 1
                error('Dimension error: Only 1 dimension currently supported')
            end
            
            % Set up spatial framework
            if isempty(model.geometry.region)
                xLims = sort(model.geometry.boundary.vertices);
                Nx = model.discretisation.Npoints;
                dx = model.discretisation.stepSize;
                xAxis = (xLims(1):dx:xLims(2))';
                op = model.discretisation.operator;
                file = model.discretisation.fileName;
                dir = model.discretisation.directory;
                [P,Q] = makeOp(op,file,Nx,dx,'directory',dir);
            else
                nRegions = length(fieldnames(model.geometry.region));
                xAxis = [];
                P = [];
                Q = [];
                for k = 1:nRegions
                    regionName = ['omega',num2str(k)];
                    xLims = sort(model.geometry.region.(regionName));
                    Nx = model.discretisation.(regionName).Npoints;
                    dx = model.discretisation.(regionName).stepSize;
                    x = (xLims(1):dx:xLims(2))';
                    xAxis = [xAxis; x];
                    op = model.discretisation.(regionName).operator;
                    file = model.discretisation.(regionName).fileName;
                    dir = model.discretisation.(regionName).directory;
                    [p,q] = makeOp(op,file,Nx,dx,'directory',dir);
                    P = blkdiag(P,p);
                    Q = blkdiag(Q,q);
                end
            end
            model.grid = xAxis;
            if ~isempty(model.solution)
                model.norm = P;
            end
            
            % Define SBP operators
            D = P\Q;
            D2 = D^2;
            
            % PDE parameters
            equationParameters = model.coefficient;
            
            % Forcing term
            if isempty(model.forcing)
                force = @(x,t) 0;
            else
                force = model.forcing;
            end
            
            % Boundary conditions
            nBC = length(fieldnames(model.boundaryCondition));
            BCpenalty = @(u,t) 0;
            for k=1:nBC
                fieldName = ['BC',num2str(k)];
                data = model.boundaryCondition.(fieldName).data;
                sigma = model.boundaryCondition.(fieldName).penalty;
                BCpoint = sscanf(model.boundaryCondition.(fieldName).node,...
                    'B%i',2);
                BCcoord = model.geometry.boundary.vertices(BCpoint);
                eVec = (xAxis==BCcoord);
                ind = find(eVec);
                
                if strcmp(model.boundaryCondition.(fieldName).BCtype,'Dirichlet')
                    BCpenalty = @(u,t) BCpenalty(u,t) + ...
                        sigma*(u - data(t)).*eVec/P(ind,ind);
                elseif strcmp(model.boundaryCondition.(fieldName).BCtype,'Neumann')
                    BCpenalty = @(u,t) BCpenalty(u,t) + ...
                        sigma*(D*u - data(t)).*eVec/P(ind,ind);
                elseif strcmp(model.boundaryCondition.(fieldName).BCtype,'Robin')
                    c = model.boundaryCondition.(fieldName).coefficients;
                    BCpenalty = @(u,t) BCpenalty(u,t) + ...
                        sigma*(c(1)*u + c(2)*D*u - data(t)).*eVec/P(ind,ind);
                end
            end
            
            % Interface conditions
            if isempty(model.geometry.interface)
                IFpenalty = @(u) 0;
            else
                nIF = length(model.geometry.interface.vertices);
                IFpenalty = @(u) 0;
                for k=1:nIF
                    fieldName = ['IF',num2str(k)];
                    sigma = model.interfaceCondition.(fieldName).penalty;
                    IFpoint = sscanf(model.interfaceCondition.(fieldName).node,...
                        'I%i',2);
                    IFcoord = model.geometry.interface.vertices(IFpoint);
                    eVec = double(sparse((xAxis==IFcoord)));
                    ind = find(eVec);
                    eVec(ind) = [sigma(1)/P(ind(1),ind(1));...
                        -sigma(2)/P(ind(2),ind(2))];
                    IFpenalty = @(u) IFpenalty(u) + ...
                        eVec*(u(ind(1)) - u(ind(2)));
                end
            end
            
            % Discretisation; u_t = RHS(t,u)
            model.RHS = @(t,u) -equationParameters(1)*u ...
                -equationParameters(2)*D*u ...
                +equationParameters(3)*D2*u ...
                +force(xAxis,t) ...
                +BCpenalty(u,t) ...
                +IFpenalty(u);
        end
    end
    
end

