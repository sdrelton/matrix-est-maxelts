function  normestm_testcode
%NORMESTM_TESTCODE  Simple tests of the codes.

% rng(3);
disp('Estimate largest entry in absolute value of A = inv(randn(100)):')
A = inv(randn(100));
[i,j,mx] = maxelt(abs(A));
[nrmestm, nrmestrow, nrmestcol, it] = normestm(A);
fprintf('Exact max = %9.2e,  ', mx)
disp(['Position = A(', num2str(i), ',', num2str(j),')'])
fprintf('Estimate  = %9.2e,  ', nrmestm)
disp(['Position = A(', num2str(nrmestrow), ',', num2str(nrmestcol),')'])

disp(' ')
disp('Also works on rectangular matrices: 100-by-50 matrix')
A = randn(100) \ rand(100,50);
[i,j,mx] = maxelt(abs(A));
[nrmestm, nrmestrow, nrmestcol, it] = normestm(A);
fprintf('Exact max = %9.2e,  ', mx)
disp(['Position = A(', num2str(i), ',', num2str(j),')'])
fprintf('Estimate  = %9.2e,  ', nrmestm)
disp(['Position = A(', num2str(nrmestrow), ',', num2str(nrmestcol),')'])

disp(' ')
disp('Now estimate largest entry (no abs value) of A = inv(randn(100)):')
A = inv(randn(100));
[i,j,mx] = maxelt(A);
opts.abs = false;
[nrmestm, nrmestrow, nrmestcol, it] = normestm(A,opts);
fprintf('Exact max = %9.2e,  ', mx)
disp(['Position = A(', num2str(i), ',', num2str(j),')'])
fprintf('Estimate  = %9.2e,  ', nrmestm)
disp(['Position = A(', num2str(nrmestrow), ',', num2str(nrmestcol),')'])

disp(' ')
disp('We can also find the largest p entries.')
if ~(exist('maxk','file')==2)
    disp(['Necessary M-file maxk is not installed, so final test will ' ...
          'not be run!'])
   disp('Type help normestm_multi for the link to maxk.')
else
   disp('Example with p = 5.')
   A = inv(randn(100));
   p = 5;
   exact = sort(abs(A(:)), 'descend');
   exact = exact(1:p);
   [nrmestm, nrmestrow, nrmestcol, it] = normestm_multi(A, p);
   disp('Compare the two lists:')
   exact'
   nrmestm
end

function [i,j,mx] = maxelt(A)
%MAXELT   Largest element of A and its indices.
[mx, ivals] = max(A);
[mx, j] = max(mx); i  = ivals(j);
