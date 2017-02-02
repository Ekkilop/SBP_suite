function [ operator ] = loadOp( strOp, fileName, directory )
%LOADOP Load SBP operator.
%   LOADOP(strOp, fileName) loads a structure containing the stencil
%   coefficients of the SBP operator with the given name (string) strOp 
%   from the file (string) fileName.
%
%   Stencil parameters are optained from the structure fields as follows:
%
%   a = operator.interior;
%   p = operator.Pboundary;
%   q = operator.Qboundary;
%
%   LOADOP(strOp, fileName, directory) operates in the directory specified
%   by (string) path.

% If no path is assigned, use current directory
if nargin == 2
    directory = cd;
end

% Check if fileName is an existing file in the path directory
if ~isFile(fileName, directory)
    error(['The file ', fileName, ' was not found in the directory ', ...
        directory])
end

% Check if strOp is an existing operator in the file
if ~isOp(strOp, fileName, directory)
    error(['The operator ', strOp, ' was not found in the file ', ...
        fileName])
end

% Include path in file name
fileName = fullfile(directory, fileName);

% Load specified operator
OP = matfile(fileName);
operatorName = matlab.lang.makeValidName(strOp);
operator = OP.(operatorName);

end