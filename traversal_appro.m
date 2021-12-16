% traversal algorithm
function [idx_A,Ps] = traversal_appro(A,k)
% to realize the CSSP by traversal method
% input: primal matrix A, maximum column k
[~,n] = size(A) ;
error_min = inf;
tic;
for c = 1:2^n-1
    if mod(c, floor(2^n/10)) == 0
        t = toc;
        clc;
        fprintf('[k=%d] %.2f(%.1f) - ', k, c/2^n, t/60);
    end
    cbin = dec2bin(c, n);
    if count(cbin, '1') == k
        idx = strfind(cbin, '1');
        [~, error_rate] = error_test(A, idx);
        if error_rate < error_min
            idx_A = idx;
            error_min = error_rate;
        end
    end
end

A = A(:,idx_A) ;
Ps = A * pinv(A) ;
