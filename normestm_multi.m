function [nrmest, nrmestrow, nrmestcol, it] = normestm_multi(A, p, opts, varargin)
% NORMESTM Estimate top-p entries of max(max(abs(A))) using matrix-vector products.
%
% The algorithm is a heavily modified version of one proposed by 
% Boyd (1974) to estimate a mixed subordinate norm,
% in our case the (1, inf)-norm to estimate max(max(abs(A))),
% using only matrix-vector products.
%
% Note: This code relies on the function maxk available at
% http://uk.mathworks.com/matlabcentral/fileexchange/23576-min-max-selection
% to find the k largest elements of an array quickly.
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
% p                   - Integer such that the top-p largest elements are
%                       estimated by the algorithm.
%
% opts (optional)     - A struct such that opts.alpha is an integer greater 
%                       than or equal to 1 and opts.abs is boolean. 
%                       The value opts.alpha controls the accuracy of the 
%                       algorithm whilst opts.abs = false will return an 
%                       estimate of max(max(A)) as opposed to 
%                       max(max(abs(A))). If opts is just an integer then 
%                       we use this value for alpha and assume 
%                       that opts.abs = true. Leaving opts empty will
%                       default to opts = 2;
%
% varargin (optional) - If A is a function handle then these additional
%                       arguments will be passed to A.
%
% ------
% Output
% ------
% nrmestm    - The estimates of max(max(abs(A))) or max(max(A)) respectively.
% nrmestmrow - The rows of the element attaining this estimate.
% nrmestmcol - The columns of the element attaining this estimate.
% it         - The number of iterations performed.
%
% -----
% Flags
% -----
% When A is a function handle the following flags will be passed to A.
%
% 'dim'      - Calling A('dim', []) returns size(A).
% 'notransp' - Calling A('notransp', x) returns A*x.
% 'transp'   - Calling A('transp', x) returns A'*x.
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
% 16th December 2015

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

try
    if nargin < 3
        alpha = 2;
        useabs = true;
    elseif isstruct(opts)
        alpha = opts.alpha;
        useabs = opts.abs;
        if useabs == 1
            useabs = true;
        end
        validateattributes(alpha, {'numeric'},...
            {'integer', 'positive', 'scalar', 'finite'})
        validateattributes(useabs, {'logical'}, {'scalar'})
    else
        useabs = true;
        alpha = opts;
        validateattributes(alpha, {'numeric'},...
            {'integer', 'positive', 'scalar', 'finite'})
    end
catch err
    error('NORMESTM:InvalidOpts',... 
        'Argument opts is invalid, please read the help for guidance.')
end
t = alpha * p;

if ismat
    [m, n] = size(A);
else
    dims = A('dim', [], varargin{:});
    m = dims(1);
    n = dims(2);
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
    X(:, 2) = transpose(b)/(3*(n/2));
end
for j = 3:t
    myrand = randi(m);
    while ismember(myrand, prev_inds)
        myrand = randi(m);
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
        [ynorms, ind] = maxk(abs(Y(:)), t);
        yvals = transpose(Y(ind));
        ynorms = ynorms';
        ind = ind';
        [ind, maxcolnum] = ind2sub(size(Y), ind);
    else
        % NOT DONE THIS YET
        [ynorms, ind] = maxk(Y(:), t);
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
            [ynormsp, indp] = maxk(abs(Yp(:)), t);
        else
            [ynormsp, indp] = maxk(Yp(:), t);
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
        [znorms, ind] = maxk(abs(Z(:)), t);
        znorms = znorms';
        ind = ind';
        [ind, ~] = ind2sub(size(Z), ind);
    else
        [znorms, ind] = maxk(Z(:), t);
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
        len = length(ind);
        if len == 0 && nrmestrow(end) ~= 0
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
% Subfunctions
%%%%%%%%%%%%%%
function [M, curindex] = ind2mat(ind, numrows, filter)
% IND2MAT Turns set of indices into a matrix of unit vectors.

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
% MATVEC Performs the matrix-vector product A*X where A may be a function.
%
% Note: We add deflation of the largest known elements throughout.
if ismat
    % Use standard matrix operations
    if transpose
        b = A'*x;
        if nrmestrow(1) ~= 0 && it > 1
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
        if nrmestrow(1) ~= 0 && it > 1
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
    len = size(x, 2);
    if transpose
        b = A('transp', x, varargin{:});
        if nrmestrow(1) ~= 0 && it > 1
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
        if nrmestrow(1) ~= 0 && it > 1
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