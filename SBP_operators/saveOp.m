function [  ] = saveOp( interior, Pboundary, Qboundary, strOp, fileName, varargin )
%SAVEOP Save SBP operator to file.
%   Given an interior stencil, (num) interior, boundary coefficients of P, 
%   (num) Pboundary, and boundary coefficients of Q, (num) Qboundary
%   SAVEOP(interior, Pboundary, Qboundary, strOp, fileName) creates a
%   structure given the name (string) strOp and saves the structure to
%   the mat-file (string) fileName.
%
%   If the file (string) fileName already exists, the new operator will be
%   appended to the file. If not, a new mat-file named (string) fileName
%   will be created.
%
%   SAVEOP(..., 'directory', directoryName) operates in the directory
%   specified by (string) directoryName.
%
%   SAVEOP(..., 'overwrite', 1) sets the flag Overwrite to True. An
%   existing operator named (string) strOp in the file (string) fileName
%   will be overwritten by a new one, specified by the input parameters.

% Check the number of input parameters
narginchk(5, 9);

% Validate input
p = inputParser;
defaultDirectory = cd;

addRequired(p,'interior',@isnumeric);
addRequired(p,'Pboundary',@isnumeric);
addRequired(p,'Qboundary',@isnumeric);
addRequired(p,'strOp',@ischar);
addRequired(p,'fileName',@ischar);
addParameter(p,'directory',defaultDirectory);
addParameter(p,'overwrite',0,...
    @(x) validateattributes(x,{'numeric'},{'integer','<=',1,'>=',0}));

parse(p, interior, Pboundary, Qboundary, strOp, fileName, varargin{:});

% Check if fileName is an existing file in the path directory
fileExists = isFile(fileName, p.Results.directory);

% If fileName exists, check if operator name is occupied
if fileExists == 1
    if isOp(strOp, fileName, p.Results.directory) && ~p.Results.overwrite
        error(['The operator name ', strOp, ' is already used in ', ...
            fileName])
    end
end

% Create field structure
interior = interior(:);
Pboundary = Pboundary(:);
Qboundary = Qboundary(:);
tempName = struct('interior',interior,'Pboundary',Pboundary,'Qboundary',Qboundary);

operatorName = matlab.lang.makeValidName(strOp);
eval([operatorName, ' = tempName;']);

% Include path in file name
fileName = fullfile(p.Results.directory, fileName);

% Save operator
if fileExists == 1
    save(fileName, eval('operatorName'), '-append')
else
    save(fileName, eval('operatorName'))
end

end

