function [varargout] = QR_appro(A,k)
% to realize the CSSP by QR method (deterministic), according to:
% [5] SUBSET SELECTION ALGORITHMS RANDOMIZED VS. DETERMINISTIC (Algrithm 1)  
% citing: [12 in article] M. Gu and S. C. Eisenstat, Efficient algorithms for computing a strong rank-revealing QR factorization,
% SIAM J. Sci. Comput., 17 (1996), pp. 848–869.
% input: primal matrix A, maximum column k

[m,n] = size(A);
% if m < n
%     error("In Method QR, m >= n is required");
% end
if k > min(m, n)
    error("In Method QR, k<=min(m,n) is required");
end
f = 1 + realmin;
[~, R, P] = qr(A);  % AP = QR
while true
    Rk = R(1:k, 1:k); Bk = R(1:k, k+1:n); Ck = R(k+1:m, k+1:n);
    Ri = inv(Rk);
    RiB = Rk\Bk;  % R_k^(-1)B_k
    gt = (RiB.^2 + sum(Ri.^2, 2) * sum(Ck.^2, 1) > f);
    if max(gt) == 0
        break
    end
    [r, c] = find(gt); i = r(1); j = c(1);
    Pi = eye(n); Pi(:, [i,j+k]) = Pi(:, [j+k,i]);
    R_tilde = R*Pi; P = P*Pi;
    [~,R] = qr(R_tilde);
end
idx_A = (1:n)*P; idx_A = idx_A(1:k);
A = A(:,idx_A) ;
Ps = A * pinv(A) ;
if nargout == 2
    varargout{1} = idx_A;
    varargout{2} = Ps;
elseif nargout == 1
    varargout{1} = P;
end