% NORMESTM_TESTCODE
rng(3);
disp('Calculate the largest entry of A = inv(randn(100))')
A = inv(randn(100));
exact = max(max(abs(A)));
[nrmestm, nrmestrow, nrmestcol, it] = normestm(A);
disp(['Exact value = ', num2str(exact)])
disp(['Estimate = ', num2str(nrmestm)])
disp(['Position = A(', num2str(nrmestrow), ',', num2str(nrmestcol),')'])

disp(' ')
disp('Also works on rectangular matrices')
A = randn(100) \ rand(100,50);
exact = max(max(abs(A)));
[nrmestm, nrmestrow, nrmestcol, it] = normestm(A);
disp(['Exact value = ', num2str(exact)])
disp(['Estimate = ', num2str(nrmestm)])
disp(['Position = A(', num2str(nrmestrow), ',', num2str(nrmestcol),')'])

disp(' ')
disp('We can also find the largest p entries. For example, p = 5.')
disp('Note: This requires the function MAXK, type help normestm_multi for the link')
A = inv(randn(100));
p = 5;
exact = sort(abs(A(:)), 'descend');
exact = exact(1:p);
[nrmestm, nrmestrow, nrmestcol, it] = normestm_multi(A, p);
disp('Compare the two lists:')
exact'
nrmestm
