classdef pdeGeometry < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimension
        boundary
        interface
        region
    end
    
    methods
        % Constructor function
        function [ geom ] = pdeGeometry( d )
            if nargin == 1
                if ~(d==1 || d==2)
                    error('Only 1 or 2 dimensions supported')
                else
                    geom.dimension = d;
                end
            end
        end
        
        % Add geometric boundary
        function addBoundary( geom, vertices, faces )
            if ~isempty(geom.boundary)
                error('Geometric error: A boundary is already defined')
            end
            
            if geom.dimension == 1
                vertices = vertices(:);
                faces = faces(:)';
                if length(vertices) ~= 2
                    error('Geometric error: Each vertex is defined by two coordinates')
                elseif length(faces) ~= 2
                    error('Geometric error: Only one face is allowed')
                else
                    geom.boundary.vertices = [vertices(1); vertices(2)];
                    geom.boundary.faces = faces;
                end
                %
            elseif geom.dimension == 2
                if size(vertices,2) ~= 2
                    error('Geometric error: Each vertex is defined by two coordinates')
                end
                geom.boundary.vertices = vertices;
                geom.boundary.faces = faces;
                %
            else
                error('Dimension error: Only 1 or 2 dimensions supported')
            end
        end
        
        % Add geometric interface
        function addInterface( geom, vertices, faces )
            
            if geom.dimension == 1
                if nargin > 2
                    error('Geometric error: No faces are allowed')
                else
                    vertices = vertices(:);
                    if isempty(geom.interface)
                        geom.interface.vertices = vertices;
                    else
                        geom.interface.vertices = [geom.interface.vertices; vertices];
                    end
                end
                %
            elseif geom.dimension == 2
                
                if size(vertices,2) ~= 2
                    error('Geometric error: Each vertex is defined by two coordinates')
                end
                if size(faces,1) ~= 2
                    error('Geometric error: Each interface must connect two points')
                end
                if isempty(geom.interface)
                    geom.interface.vertices = vertices;
                    geom.interface.faces = faces;
                else
                    geom.interface.vertices = [geom.interface.vertices; vertices];
                    geom.interface.faces = [geom.interface.faces; faces];
                end
                %
            else
                error('Dimension error: Only 1 or 2 dimensions supported')
            end
            
        end
        
        % Get boundary vertex coordinates
        function [ ind ] = getBVertex( geom, vertex )
            if length(vertex) ~= 2 || ~ischar(vertex) || vertex(1) ~= 'B'
                error('Input error: Expected string input of form Bn')
            end
            ind = geom.boundary.vertices(str2double(vertex(2)),:);
        end
        
        % Get interface vertex coordinates
        function [ ind ] = getIVertex( geom, vertex )
            if length(vertex) ~= 2 || ~ischar(vertex) || vertex(1) ~= 'I'
                error('Input error: Expected string input of form In')
            end
            ind = geom.interface.vertices(str2double(vertex(2)),:);
        end
        
        % Define geometric region
        function defineRegion( geom, nodes )
            if ~iscell(nodes)
                error('Input error: Cell input expected')
            end
            
            nNodes = length(nodes);
            
            if geom.dimension == 1 && nNodes ~= 2
                error('Geometric error: Two nodes expected')
            end
            if geom.dimension == 2 && nNodes ~= 4
                error('Geometric error: Four nodes expected')
            end
            
            if isempty(geom.region)
                regionName = 'omega1';
            else
                nNames = length(fieldnames(geom.region));
                regionName = ['omega', num2str(nNames+1)];
            end
            geom.region.(regionName) = [];
            for k = 1:nNodes
                if nodes{k}(1) == 'B'
                   geom.region.(regionName) = [geom.region.(regionName);...
                       getBVertex(geom,nodes{k})];
                elseif nodes{k}(1) == 'I'
                    geom.region.(regionName) = [geom.region.(regionName);...
                        getIVertex(geom,nodes{k})];
                else
                    error('Input error: Inappropriate node index')
                end
            end

        end
        
        % Plot geometry
        function plotGeometry( geom )
            
            % Set plot ranges
            xMin = min(geom.boundary.vertices(:,1));
            xMax = max(geom.boundary.vertices(:,1));
            xRange = xMax - xMin;
            xMin = xMin - 0.1*xRange;
            xMax = xMax + 0.1*xRange;
            xRange = xMax - xMin;
            
            % Draw boundaries
            if geom.dimension == 1
                iPoints = geom.boundary.vertices;
                line([iPoints(1) iPoints(2)],[0 0],'Color','black',...
                    'LineWidth',2,'Marker','o','MarkerSize',8,...
                    'MarkerFaceColor','black')
                yMin = -1;
                yMax = 1;
                yRange = yMax - yMin;                
                %
            else
            patch('Vertices',geom.boundary.vertices,'Faces',geom.boundary.faces,...
                'LineWidth',2,'FaceColor','none','Marker','o','MarkerSize',8,...
                'MarkerFaceColor','black');
            yMin = min(geom.boundary.vertices(:,2));
            yMax = max(geom.boundary.vertices(:,2));
            yRange = yMax - yMin;
            yMin = yMin - 0.1*yRange;
            yMax = yMax + 0.1*yRange;
            yRange = yMax - yMin;
            end
            axis([xMin xMax yMin yMax])
            
            % Annotate boundary points
            axisPos = get(gca,'position');
            axisXMin = axisPos(1);
            axisXMax = axisPos(1) + axisPos(3);
            axisYMin = axisPos(2);
            axisYMax = axisPos(2) + axisPos(4);
            xRangeScaled = axisXMax - axisXMin;
            yRangeScaled = axisYMax - axisYMin;
            for k = 1:size(geom.boundary.vertices,1)
                xPos = geom.boundary.vertices(k,1);
                if geom.dimension == 1
                    yPos = 0;
                else
                    yPos = geom.boundary.vertices(k,2);
                end
                xPosScaled = (xPos - xMin)*xRangeScaled/xRange + axisXMin;
                yPosScaled = (yPos - yMin)*yRangeScaled/yRange + axisYMin;
                dim = [xPosScaled yPosScaled 0 0];
                str = ['$B',num2str(k),'$'];
                t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
                t.Interpreter = 'latex';
                t.FontSize = 14;
                t.EdgeColor = 'none';
            end
            
            % Draw interfaces
            if ~isempty(geom.interface)
                if geom.dimension == 1
                    iPoints = geom.interface.vertices;
                    for k = 1:length(iPoints)
                        line([iPoints(k) iPoints(k)],[-0.1 0.1],...
                            'Color','black')
                    end
                    %
                else
                    patch('Vertices',geom.interface.vertices,...
                        'Faces',geom.interface.faces,...
                        'FaceColor','none','LineStyle','--',...
                        'EdgeColor','red',...
                        'Marker','o','MarkerSize',8,...
                        'MarkerFaceColor','red',...
                        'MarkerEdgeColor','black')
                end
                
                % Annotate interface points
                for k = 1:size(geom.interface.vertices,1)
                    if geom.dimension == 1
                        xPos = geom.interface.vertices(k);
                        yPos = 0.1;
                    else
                        xPos = geom.interface.vertices(k,1);
                        yPos = geom.interface.vertices(k,2);
                    end
                    xPosScaled = (xPos - xMin)*xRangeScaled/xRange + axisXMin;
                    yPosScaled = (yPos - yMin)*yRangeScaled/yRange + axisYMin;
                    yPosScaled = yPosScaled + 0.06*yRangeScaled;
                    dim = [xPosScaled yPosScaled 0 0];
                    str = ['$I',num2str(k),'$'];
                    t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
                    t.Interpreter = 'latex';
                    t.FontSize = 14;
                    t.EdgeColor = 'none';
                    t.Color = 'red';
                end
                
            end
            
            % Draw regions
            % 2D
            if ~isempty(geom.region) && geom.dimension == 2
                nRegions = length(fieldnames(geom.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    nodes= geom.region.(fieldName);
                    nNodes = size(nodes,1);
                    shade = k/nRegions;
                    p = patch('Vertices',nodes,'Faces',1:nNodes,...
                        'EdgeColor','None',...
                        'FaceColor',[0.5 0.5 0.5],...
                        'FaceAlpha',shade);
                    uistack(p,'bottom')
                    
                    % Annotate regions
                    xPos = (max(nodes(:,1)) - min(nodes(:,1)))/2 + min(nodes(:,1));
                    yPos = (max(nodes(:,2)) - min(nodes(:,2)))/2 + min(nodes(:,2));
                    xPosScaled = (xPos - xMin)*xRangeScaled/xRange + axisXMin;
                    yPosScaled = (yPos - yMin)*yRangeScaled/yRange + axisYMin;
                    xPosScaled = xPosScaled - 0.025*xRangeScaled;
                    yPosScaled = yPosScaled + 0.025*yRangeScaled;
                    dim = [xPosScaled yPosScaled 0 0];
                    str = ['$\Omega_',num2str(k),'$'];
                    t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
                    t.Interpreter = 'latex';
                    t.FontSize = 16;
                    t.EdgeColor = 'none';
                    t.Color = 'black';
                end
                % 1D
            elseif ~isempty(geom.region) && geom.dimension == 1
                nRegions = length(fieldnames(geom.region));
                for k = 1:nRegions
                    fieldName = ['omega',num2str(k)];
                    nodes= geom.region.(fieldName);
                    
                    % Annotate regions
                    xPos = (max(nodes) - min(nodes))/2 + min(nodes);
                    yPos = -0.1;
                    xPosScaled = (xPos - xMin)*xRangeScaled/xRange + axisXMin;
                    yPosScaled = (yPos - yMin)*yRangeScaled/yRange + axisYMin;
                    xPosScaled = xPosScaled - 0.025*xRangeScaled;
                    dim = [xPosScaled yPosScaled 0 0];
                    str = ['$\Omega_',num2str(k),'$'];
                    t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
                    t.Interpreter = 'latex';
                    t.FontSize = 16;
                    t.EdgeColor = 'none';
                    t.Color = 'black';
                end
            end
            
        end
    end
    
end

