function [elements, index] = maxk_default(v, k)
%maxk_default Find the k largest elements of the vector v.
% This function calculates the largest k elements of the input vector v.
% It is designed primarily for use with normestm and normestm_multi.
%
% -----
% Input
% -----
% v - Vector of real number entries.
% k - Integer stating how many of the largest elements are to be found.
%
% ------
% Output
% ------
% elements - List of the largest k elements of v.
% index    - Indices where the largest k elements in v can be found.
%
% ------
% Author
% ------
% Samuel D. Relton
% 19th December 2015

mycopy = v;
index = zeros(k, 1);

for i = 1:k
    [~, idx] = max(mycopy);
    mycopy(idx) = -inf;
    index(i) = idx;
end
elements = v(index);

end