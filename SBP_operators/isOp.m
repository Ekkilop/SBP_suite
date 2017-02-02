function [ opExists ] = isOp( strOp, fileName, directory )
%ISOP Check for operator in mat-file.
%   ISOP(strOp, fileName) returns 1 if the operator with the given name 
%   (string) strOp exists in the file (string) fileName and 0 otherwise.
%
%   ISOP(strOp, fileName, directory) operates in the directory specified by
%   (string) directory.

% If no path is assigned, use current directory
if nargin == 2
    directory = cd;
end

% If fileName exists, search for operator
if isFile(fileName, directory)
    operatorName = matlab.lang.makeValidName(strOp);
    operatorList = who('-file', fileName);
    opExists = ismember(operatorName, operatorList);
    
    % If fileName does not exist, return 0
else
    warning(['The file ', fileName, ' does not exist in the current' ...
        ' directory. isOp.m returns 0 automatically.'])
    opExists = 0;
end
end