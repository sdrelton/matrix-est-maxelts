function [nrmest, nrmestrow, nrmestcol, it] = normestm_multi(A, p, opts, varargin)
%normestm_multi Estimate largest p entries of A or abs(A).
% This function estimates the largest p entries of A or abs(A).
% The estimate is computed using only matrix-vector products with A and A'.
% The general syntax of a function call is
%    [nrmest, nrmestrow, nrmestcol, it] = normestm_multi(A, p, opts, varargin)
%
% Note: This code relies on the function maxk_default.
% A more optimized version is available at
% http://uk.mathworks.com/matlabcentral/fileexchange/23576-min-max-selection
% to find the k largest elements of an array quickly.
% It is not used directly here due to compilation issues on some platforms.
% To use maxk, replace all occurences of maxk_default with maxk.
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
% p                   - Integer such that the largest p elements are
%                       estimated by the algorithm.
%
% opts (optional)     - A structure with two fields.
%                       opts.alpha: an integer greater than or equal to 1, 
%                          which control the number of columns in the matrix 
%                          iterates.  Larger values generally give more 
%                          accurate estimates.
%                        opts.abs: a logical.  If true, max(max(abs(A)))
%                          is estimated, otherwise max(max(A)).
%                        The defaults are opts.alpha = 2, opts.abs = true.    
%
% varargin (optional) - If A is a function handle then these additional
%                       arguments will be passed to A.
%
% ------
% Output
% ------
% nrmestm    - The estimates of max(max(abs(A))) or max(max(A)),
%              depending on opts.abs.
% nrmestmrow - The rows of the element attaining this estimate.
% nrmestmcol - The columns of the element attaining this estimate.
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
% Corresponding algorithm from the paper: Algorithm 5.2
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
    error('NORMESTM_MULTI:NoInput', 'No inputs to normestm_multi.')
end
ismat = isnumeric(A);
if ismat
    validateattributes(A, {'numeric'}, {'2d'})
else
    validateattributes(A, {'function_handle'}, {})
end

if nargin < 2
    error('NORMESTM_MUTLI:NotEnoughInputs', 'Need at least 2 inputs.')
end

% Defaults.
alpha = 2;
useabs = true;

if nargin >= 3
   try
       if isstruct(opts)
           if isfield(opts,'alpha')
              alpha = opts.alpha;
           end   
            if isfield(opts,'abs')
               useabs = logical(opts.abs);
           end
       else
           alpha = opts;
       end
        validateattributes(alpha, {'numeric'},...
            {'integer', 'positive', 'scalar', 'finite'})
   catch err
       error('NORMESTM:InvalidOpts',... 
           'Argument opts is invalid, please read the help for guidance.')
   end
end
t = alpha * p;

if ismat
    [m, n] = size(A);
else 
    dims = A('dim', [], varargin{:});
    m = dims(1); n = dims(2);
end

% Initialise the algorithm
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
    clear X
    % Find largest elements of Y
    if useabs
        [ynorms, ind] = maxk_default(abs(Y(:)), t);
        yvals = transpose(Y(ind));
        ynorms = ynorms';
        ind = ind';
        [ind, maxcolnum] = ind2sub(size(Y), ind);
    else
        [ynorms, ind] = maxk_default(Y(:), t);
        yvals = transpose(Y(ind));
        ynorms = ynorms';
        ind = ind';
        [ind, maxcolnum] = ind2sub(size(Y), ind);
    end
    maxynorm = ynorms;
    if it == 1 && t > 2
        % Store best elements from first iteration
        Yp = Y(:, 3:end);
        if useabs
            [ynormsp, indp] = maxk_default(abs(Yp(:)), t);
        else
            [ynormsp, indp] = maxk_default(Yp(:), t);
        end
        actualvals = transpose(Yp(indp(1:p)));
        ynormsp = ynormsp';
        indp = indp';
        [indp, maxcolnump] = ind2sub(size(Y), indp);
        nrmestcol = prev_inds(maxcolnump(1:p));
        nrmestrow = indp(1:p);
        nrmest = ynormsp(1:p);
        clear Yp
    elseif it > 1
        maxcol = curcols(maxcolnum);
        maxrow = ind;
        if any((maxynorm > nrmest(end)) & ~(ismember(maxrow, nrmestrow) & ismember(maxrow, nrmestcol)))
            % Update largest elements found so far
            nrmest = [nrmest, maxynorm];
            nrmestrow = [nrmestrow, maxrow];
            nrmestcol = [nrmestcol, maxcol];
            actualvals = [actualvals, yvals];
            % Now double check we haven't got the same element twice.
            U = [nrmest', nrmestrow', nrmestcol', actualvals'];
            [~, uniquerows] = unique(U(:,[2,3]) ,'rows');
            uniquerows = uniquerows(end:-1:1);
            U = U(uniquerows, :);
            [~, order] = sort(U(:,1), 'descend');
            U = U(order, :);
            % Now possible that U has less than p rows.
            nrmest = U(1:min(p, size(U, 1)), 1)';
            nrmestrow = U(1:min(p, size(U, 1)), 2)';
            nrmestcol = U(1:min(p, size(U, 1)), 3)';
            actualvals = U(1:min(p, size(U, 1)), 4)';
        elseif nrmestrow(end) ~= 0
            % No better entry found so quit
            break
        end
    end
    clear Y
    
    % Compute Z = A'*dual(Y)
    Z = matvec(A, ind2mat(ind, m, false), ismat, true);
    if useabs
        [znorms, ind] = maxk_default(abs(Z(:)), t);
        znorms = znorms';
        ind = ind';
        [ind, ~] = ind2sub(size(Z), ind);
    else
        [znorms, ind] = maxk_default(Z(:), t);
        znorms = znorms';
        ind = ind';
        [ind, ~] = ind2sub(size(Z), ind);
    end
    clear Z
    if it > 1
        maxznorm = max(znorms);
        if all(maxznorm <= maxynorm) && nrmestrow(end) ~= 0
            % Converged so exit
            break
        end
        ind = setdiff(ind, prev_inds);
        if isempty(ind) && nrmestrow(end) ~= 0
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

%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%
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
% Deflation of the largest known elements is used.
if ismat
    % Use standard matrix operations
    if transpose
        b = A'*x;
        if it > 1
            for k = 1:length(nrmestrow)
                if nrmestrow(k) == 0 
                    continue
                else
                    Aji = conj(actualvals(k));
                    b(nrmestcol(k), :) = b(nrmestcol(k), :) - Aji * x(nrmestrow(k), :);
                end
            end
        end
    else
        b = A*x;
        if it > 1
            for k = 1:length(nrmestrow)
                if nrmestrow(k) == 0 
                    continue
                else
                    Aij = actualvals(k);
                    b(nrmestrow(k), :) = b(nrmestrow(k), :) - Aij * x(nrmestcol(k), :);
                end
            end
        end
    return
    end
else
    % A is a function, use API described in help text.
    if transpose
        b = A('transp', x, varargin{:});
        if it > 1
            for k = 1:length(nrmestrow)
                if nrmestrow(k) == 0 
                    continue
                else                
                    Aji = conj(actualvals(k));
                    b(nrmestcol(k), :) = b(nrmestcol(k), :) - Aji*x(nrmestrow(k), :);
                end
            end
        end
    else
        b = A('notransp', x, varargin{:});
        if it > 1
            for k = 1:length(nrmestrow)
                if nrmestrow(k) == 0
                    continue
                else                
                    Aij = actualvals(k);
                    b(nrmestrow(k), :) = b(nrmestrow(k), :) - Aij*x(nrmestcol(k), :);
                end
            end
        end
    end
    return       
end
end % End of matvec

end