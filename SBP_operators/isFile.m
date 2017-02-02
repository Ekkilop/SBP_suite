function [ fileExists ] = isFile( fileName, directory )
%ISFILE Cheack for file in directory
%   ISFILE(fileName) returns 1 if a file with the given name (string)
%   fileName exists in the current directory and 0 otherwise.
%
%   ISFILE(fileName, directory) operates in the directory specified by
%   (string) directory.

% If path is not given, search in current directory
if nargin == 1
    directory = cd;
end

% Search for file
if exist(directory,'dir') == 7
    if exist(fullfile(directory,fileName),'file') == 2
        fileExists = 1;
    else
        fileExists = 0;
    end
else
    error(['The directory ' directory ' was not found on the search path.'])
end


end