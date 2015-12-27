# matrix-est-maxelts
MATLAB code to estimate the largest elements of a matrix using only matrix-vector products.

## Introduction
The functions in this repository are designed to estimate
the largest elements of A or abs(A), 
where the m-by-n matrix A is not known explicitly.  For example we may have A =
B*C, A = expm(B), or A = inv(B), where forming A explicitly is impractical
(maybe even impossible if B and C are large and sparse).  However in each
case forming matrix-vector products with A and its conjugate transpose may
nevertheless be possible and even relatively inexpensive.

The function `normest` estimates 
(in MATLAB notation)
- `max(max(abs(A)))` (the default), or
- `max(max(A))`.

The quantity `max(max(abs(A)))` can be expressed as the mixed subordinate
(1,inf)-norm of A.  Underlying our algorithms is an algorithm of Boyd
(1974) and Tao (1975) for estimating mixed subordinate norms.

The function `normestm_multi` estimates the largest p elements of A or abs(A),
where p is an input argument.
The function depends upon maxk_default.
A more optimized version of this function is available,
authored by Bruno Luong.
This free software can be downloaded from the
[MATLAB File
Exchange](http://uk.mathworks.com/matlabcentral/fileexchange/23576-min-max-selection)
but it uses MEX files and can be difficult to install and get working.
To make use of this code, should you be able to install and run it successfully,
replace all occurrences of maxk_default with maxk in normestm_multi.m

Full details of the algorithms,
along with thorough numerical experiments
investigating their performance, can be found in the (open access) paper

N. J. Higham and S. D. Relton, "[Estimating the Largest Elements of a Matrix](http://eprints.ma.man.ac.uk/2424)", MIMS EPrint
2015.116, Manchester Institute for Mathematical Sciences, The University of
Manchester, UK, December 2015.

To check that the code is functioning properly you can run
[the testcode](normestm_testcode.m) in MATLAB.

## Details
Our algorithms can be applied to matrices known explicitly,
which is usually less efficient than just calling max(max(abs(A)))
directly, or to matrices known only implicitly. Here is an example of the
explicit case.  The implicit case is discussed below.

```matlab
A = inv(randn(500));
opts.t = 3;
opts.abs = true;
[nrmest, nrmrow, nrmcol, iters] = normestm(A, opts);

clear all
A = inv(randn(500));
opts.alpha = 5;
opts.abs = true;
p = 5;
[nrmests, nrmrows, nrmcols, iters] = normestm_multi(A, p, opts);
```
##### Normest inputs and outputs
The inputs to `normestm` are:
* `A`:     Matrix with real/complex entries or a function handle (see the advanced examples below). Can be rectangular.
* `opts`:  Structure where `opts.t` (`opts.alpha` for `normestm_multi`) is an integer
       controlling the accuracy and `opts.abs` is a logical variable that
       determines whether to seek largest elements of `max(max(abs(A)))` or `max(max(A))`. When `opts` is
       an integer, instead of a structure, `opts.abs` is set to true.

In addition `normestm_multi` has one extra input:
* `p` - An integer denoting how many of the largest elements are required.

The outputs of `normestm` and `normestm_multi` are:
* `nrmest:`     The estimate(s) of the largest p elements.
* `nrmestrow`:  The rows where the estimates appear.
* `nrmestcol`:  The columns where the estimates appear.
* `iter`:       The number of iterations required.

##### Using implicitly defined matrices
The main power of this algorithm is that it can be applied to matrices
where only matrix-vector products can be computed.  To make use of this
functionality you will need to write a short wrapper function to perform
the matrix-vector products in the following form.

```matlab
function b = mywrapper(flag, x)
%mywrapper Computes matrix-vector products for use with normestm and normestm_multi.

switch flag
    case 'dim'
        % Compute size of the implicit matrix (number of rows and columns).
        b = ...;
    case 'notransp'
        % Compute A*x
        b = ...
    case 'transp'
        % Compute A'*x
        b = ...
    otherwise
        error('Undefined flag');
end
```

A fully functioning wrapper for computing the matrix exponential (using [expmv](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector))
is included in the repository as [`expmv_wrapper.m`.](expmv_wrapper.m)
We can find the p = 10 largest entries of the matrix exponential as follows
(the matrix is taken from the University of Florida Sparse Matrix Collection).
Here we use the function UFget, which can be downloaded from
[MATLAB File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/11896-ufget--matlab-interface-to-the-uf-sparse-matrix-collection#comments).

```matlab
p = 10;
opts.alpha = 3;
opts.abs = true;
A = UFget('SNAP/ca-AstroPh');
A = A.A;
[nrmests, nrmows, nrmcols] = normestm_multi(@expmv_wrapper, p, opts, A);
```
