% select subset of columns by random selecting
function [idx_A,Ps] = random_appro(A,k)
% to realize the CSSP by random choosing method
% input: primal matrix A, maximum column k
[~,n] = size(A) ;
idx_A = randperm(n, k);
sort(idx_A);
A = A(:,idx_A) ;
Ps = A * pinv(A) ;
