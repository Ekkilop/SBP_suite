classdef pdeDiscretisation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimension
        Npoints
        operator
        stepSize
    end
    
    properties (Hidden)
        fileName
        directory
    end
    
    methods
        % Constructor function
        function [ disc ] = pdeDiscretisation( d )
            if nargin == 1
                if ~(d==1 || d==2)
                    error('Only 1 or 2 dimensions supported')
                else
                    disc.dimension = d;
                end
            end
        end
        
        % Set resolution
        function setRes( disc, nPoints )
            if disc.dimension == 1 && numel(nPoints) ~= 1
                error('Input error: Scalar input expected')
            elseif disc.dimension == 2 && numel(nPoints) ~= 2
                error('Input error: Vector with two elements required')
            end
            disc.Npoints = nPoints(:);
        end
        
        % Set operator
        function setOp( disc, op, file, dir )
            if ~(ischar(op) && ischar(file))
                error('Input error: String input expected')
            end
            disc.operator = op;
            disc.directory = dir;
            if nargin == 4
                disc.fileName = file;
            else
                disc.directory = cd;
            end
        end
        
        % Calculate step
        function calcStep( disc, region )
            region = sort(region);
            span = region(end) - region(1);
            disc.stepSize = span/disc.Npoints;
        end
    end
end

