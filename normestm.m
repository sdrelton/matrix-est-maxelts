function [nrmest, nrmestrow, nrmestcol, it] = normestm(A, opts, varargin)
%normestm Estimate largest entry of A or abs(A) using matrix-vector products.
% This function estimates the largest entry of A or abs(A).
% The latter quantity is equal to mixed subordinate (1, inf)-norm of A.
% The estimate is computed using only matrix-vector products with A and A'.
% The general syntax of a function call is
%   function [nrmest, nrmestrow, nrmestcol, it] = normestm(A, opts, varargin)
%
% -----
% Input
% -----
% A                   - Matrix with real or complex entries, or a 
%                       function handle such that A('notransp', x) 
%                       computes A*x and A('transp', x) computes A'*x.
%                       Please see below for additional requirements on the 
%                       call A(x, flag).
%
% opts (optional)     - A structure with two fields.
%                       opts.t: an integer greater than or equal to 1, 
%                          which is the number of columns in the matrix 
%                          iterates.  Larger values generally give more 
%                          accurate estimates.
%                        opts.abs: a logical.  If true, max(max(abs(A)))
%                          is estimated, otherwise max(max(A)).
%                        The defaults are opts.t = 2, opts.abs = true.    
%
% varargin (optional) - If A is a function handle then these additional
%                       arguments will be passed to A.
%
% ------
% Output
% ------
% nrmestm    - An estimate of max(max(abs(A))) or max(max(A)), depending
%              on opts.abs.
% nrmestmrow - The row of the element attaining this estimate.
% nrmestmcol - The column of the element attaining this estimate.
% it         - The number of iterations performed.
%
% -----
% Flags
% -----
% When A is a function handle the following flags will be passed to A.
%
% 'dim'      - Calling A('dim', []) should return size(A).
% 'notransp' - Calling A('notransp', x) should return A*x.
% 'transp'   - Calling A('transp', x) should return A'*x.
%
% ---------
% Reference
% ---------
% Nicholas J. Higham and Samuel D. Relton,
% Estimating the Largest Elements of a Matrix,
% MIMS EPrint ....., 2015.
%
% Corresponding algorithm from the paper: Algorithm 4.1
%
% -------
% Authors
% -------
% Nicholas J. Higham and Samuel D. Relton
% 18th December 2015

%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%
% Check all input arguments
if nargin == 0
    error('NORMESTM:NoInput', 'No inputs to normestm.')
end
ismat = isnumeric(A);
if ismat
    validateattributes(A, {'numeric'}, {'2d'})
else
    validateattributes(A, {'function_handle'}, {})
end

% Defaults.
t = 2;
useabs = true;

if nargin >= 2
   try
       if isstruct(opts)
           if isfield(opts,'t')
              t = opts.t;
           end   
            if isfield(opts,'abs')
               useabs = logical(opts.abs);
           end
       else
           t = opts;
       end
        validateattributes(t, {'numeric'},...
            {'integer', 'positive', 'scalar', 'finite'})
   catch err
       error('NORMESTM:InvalidOpts',... 
           'Argument opts is invalid, please read the help for guidance.')
   end
end

if ismat
    [m, n] = size(A);
else
    dims = A('dim', [], varargin{:});
    m = dims(1); n = dims(2);
end

% Initialize the algorithm.
maxiter = 20;
prev_inds = [];
nrmest = 0;
nrmestrow = 0;
nrmestcol = 0;

X = zeros(n, t);
X(:, 1) = ones(n,1)/n;
if t >  1
    b = 1 + ((1:n) - 1)/(n-1);
    b(2:2:end) = -b(2:2:end);
    X(:, 2) = transpose(b)/(3*n/2);
end
for j = 3:t
    myrand = randi(n);
    while ismember(myrand, prev_inds)
        myrand = randi(n);
    end
    X(myrand, j) = 1;
    prev_inds(end+1) = myrand;
end

%%%%%%%%%%%
% Main loop
%%%%%%%%%%%
for it = 1:maxiter
    % Compute Y = A*X
    Y = matvec(A, X, ismat, false);
    % Find largest elements of Y
    if useabs
        [ynorms, ind] = max(abs(Y));
    else
        [ynorms, ind] = max(Y);
    end
    [maxynorm, maxcolnum] = max(ynorms);
    if it == 1 && t > 2
        % Store best element from first iteration
       [nrmest, nrmestcol] = max(ynorms(3:end));
       nrmestrow = prev_inds(nrmestcol);
    elseif it > 1
        if maxynorm > nrmest
            % Update largest element found so far
            maxcol = curcols(maxcolnum);
            maxrow = ind(maxcolnum);
            nrmest = maxynorm;
            nrmestrow = maxrow;
            nrmestcol = maxcol;
        else
            % No better entry found so quit
            break
        end
    end
    
    % Compute Z = A'*dual(Y)
    Z = matvec(A, ind2mat(ind, m, false), ismat, true);
    if useabs
        [znorms, ind] = max(abs(Z));
    else
        [znorms, ind] = max(Z);
    end
    if it > 1
        maxznorm = max(znorms);
        if maxznorm <= maxynorm
            % Converged so exit
            break
        end
        ind = setdiff(ind, prev_inds);
        if length(ind) == 0
            break
        end
    end % End of convergence test
    % Get unit vectors for next iterations
    [X, curcols] = ind2mat(ind, n, true);
    if curcols == -1
        % Not enough columns left to do next iteration
        break;
    end
    % Update which unit vectors we have used
    prev_inds = [prev_inds, ind];
end % End of main loop

%%%%%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%%%%%%

function [M, curindex] = ind2mat(ind, numrows, filter)
%ind2mat Turns set of indices into a matrix of unit vectors.

if filter
    % Find unused indices
    ind = setdiff(ind, prev_inds);
end
len = length(ind);
curindex = zeros(1, t);
curindex(1:len) = ind;

M = zeros(numrows, t);
M(sub2ind([numrows, len], ind, 1:len)) = 1;
if filter
    rem_cols = t - len;
    if rem_cols > 0
        remindex = setdiff(1:n, prev_inds);
        if length(remindex) < rem_cols
            % Not enough columns for next iteration, quit now
            curindex = -1;
            return;
        end
        for k = 1:rem_cols
            index = randi(length(remindex));
            curindex(len+k) = remindex(index);
            M(remindex(index), len+k) = 1;
            remindex = [remindex(1:index-1), remindex(index+1:end)];
        end
    
    end
end
end % End of ind2mast

function b = matvec(A, x, ismat, transpose)
%matvec Performs the matrix-vector product A*X where A may be a function.
if ismat
    % Use standard matrix operations
    if transpose
        b = A'*x;
    else
        b = A*x;
    end
    return
else
    % A is a function, use API described in help text.
    if transpose
        b = A('transp', x, varargin{:});
    else
        b = A('notransp', x, varargin{:});
    end
    return       
end
end

end