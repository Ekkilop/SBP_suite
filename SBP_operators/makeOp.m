function [ P, Q ] = makeOp( strOp, fileName, N, dx, varargin )
%MAKEOP creates an SBP operator.
%   [P, Q] = MAKEOP(strOp, fileName, N, dx) creates an SBP operator
%   with the stencil data provided by the operator (string) strOp,
%   loaded from the file (string) fileName, and
%   the problem specific parameters (num) N and (num) dx.
%
%   MAKEOP assumes that the full operator should be created.
%   If the periodic operator is sought, use 
%   MAKEOP(..., 'form', 'periodic').
%
%   MAKEOP(..., 'directory', directoryName) operates in the directory
%   specified by (string) directoryName.

% Check the number of input parameters
narginchk(4, 9);

% Validate input
par = inputParser;
defaultDirectory = cd;
defaultForm = 'full';
expectedForms = {'full', 'periodic'};

addRequired(par,'strOp',@ischar);
addRequired(par,'fileName',@ischar);
addRequired(par,'N',@isnumeric);
addRequired(par,'dx',@isnumeric);
addParameter(par,'directory',defaultDirectory);
addParameter(par,'form',defaultForm,...
    @(x) any(validatestring(x,expectedForms)));

parse(par, strOp, fileName, N, dx, varargin{:});

% Check if the operator exists
if ~isOp(strOp, fileName, par.Results.directory)
    error('The sought operator does not exist.')
end

% Load specified operator
operator = loadOp(strOp, fileName, par.Results.directory);

% Extract stencil parameters
a = operator.interior;
p = operator.Pboundary;
q = operator.Qboundary;

% Semi-bandwidth and boundary dimension
width = length(a);
r = length(p);

% Periodic operator
% {
if strcmp(par.Results.form, 'periodic')
    P = dx*speye(N);

    diagQ = repmat([-fliplr(a'), 0, a'], N, 1);
    Q = spdiags(diagQ,-width:width,N,N) + ...
        spdiags(diagQ(:,1:width),(N-width):(N-1),N,N) + ...
        spdiags(diagQ(:,width+2:end),-(N-1):-(N-width),N,N);
%}
% Full operator    
elseif strcmp(par.Results.form, 'full')
    P = sparse(dx*diag([p', ones(1,N-(2*r-1)), fliplr(p')]));
    
    diagQ = repmat([-fliplr(a'), 0, a'], N+1, 1);
    Q = spdiags(diagQ,-width:width,N+1,N+1);
    
    Qbound = tril(ones(r),-1);
    Qbound(Qbound == 1) = -q;
    Qbound = Qbound - Qbound';
    Qbound(1,1) = -1/2;
    
    Q(1:r,1:r) = Qbound;
    Q((N-r+2):end,(N-r+2):end) = -rot90(Qbound,2);
end

end

