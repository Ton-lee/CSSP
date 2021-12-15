function [idx_A,Ps] = randomized_appro(A,k)
% to realize the CSSP by randomized method, according to:
% [5] SUBSET SELECTION ALGORITHMS RANDOMIZED VS. DETERMINISTIC (Algrithm 1,2,3)
% citing [3, 1 in article] An Improved Approximation Algorithm for the Column Subset Selection Problem
% input: primal matrix A, maximum column k
[~,n] = size(A);
[~, S, V] = svd(A,'econ');
Vk = V(:, 1:k);
S_ = S(k+1:end, k+1:end); V_ = V(:, k+1:end);
p = sum(Vk.^2, 2)'/(2*k) + sum((S_*V_').^2, 1)/sum(sum((S_*V_').^2, 1))/2;
error_min = inf;
% 40 repeats
for repeat = 1:40
    c = k;
    % choose c
    S_W = zeros(1, k);
    while S_W(k) < 1/2 && c <= n
        c = 2*c;
        % randomized stage
        D = diag(1./sqrt(min(1, c*p)));
        P = eye(n); W = Vk'*D;
        c_tilde = 0;
        for i = 1:n
            q = rand();
            if q < c*p(i)
                c_tilde = c_tilde + 1;
                Pi = eye(n); Pi(:, [i,c_tilde]) = Pi(:, [c_tilde,i]);
                W = W*Pi; P = P*Pi;
            end
        end
        WR = W(:, 1:c_tilde);
        S_W = svd(WR);
        if length(S_W) < k  % in case that there are not enough column selected
            S_W = zeros(1, k);
        end
    end
    % deterministic stage
    if length(WR) <= k  % unluckily we get too little columns
        Pi = eye(c_tilde);
    else
        Pi = QR_appro(WR, k);
    end
    P = P * [[Pi, zeros(c_tilde, n-c_tilde)];[zeros(n-c_tilde, c_tilde), eye(n-c_tilde)]];
    idx = (1:n)*P; idx = idx(1:k);
    [~, error_rate] = error_test(A, idx);
    if error_rate < error_min
        idx_A = idx;
        error_min = error_rate;
    end
end
% result
A = A(:,idx_A);
Ps = A * pinv(A);