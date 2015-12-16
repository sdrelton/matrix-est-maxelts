function b = expmv_wrapper(flag, x, A)
% EXPMV_WRAPPER Wrapper for use with normestm to get max(max(abs(expm(A))).

% The code for expmv can be found at:
% http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector

switch flag
    case 'dim'
        b = size(A);
    case 'notransp'
        b = expmv(1.0, A, x);
    case 'transp'
        b = expmv(1.0, A', x);
    otherwise
        error('Undefined flag');
end
    

end