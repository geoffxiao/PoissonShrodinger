% Null space = columns of N
N = [1 2 3 4; 5 6 7 8]'
N = [-2 -2 -2 -1; 3 3 3 1]'

[r c] = size(N);
% zero rows = r - c
% c = dim of null space
% r = dim of R^n vector space

A = [rref(N'); zeros(r - c, r)]