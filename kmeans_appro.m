function [idx_A,Ps] = kmeans_appro(A,k,choice)
% to realize the CSSP by kmeans method
% input: primal matrix A, maximum column k, choice
% choice == 1 : k clusters + nearest neighbor
% choice == 2 : delete columns one by one
[~,n] = size(A) ;
A_save = A ; % copy A
% normalization by 2-norm
for a = 1:n
    A(:,a) = A(:,a) / norm(A(:,a),2) ;
end
A = A.' ;
if choice == 1
    idx_A = zeros(1,k) ;
    % k clusters and the distances
    [clus_idx,~,~,D1] = kmeans(A,k) ;
    D = zeros(1,n) ;
    for temp_d = 1:n
        D(temp_d) = D1(temp_d,clus_idx(temp_d)) ;
    end
    % find k columns nearest to the center
    for a = 1:k
        idx_clus = find(clus_idx == a) ;
        [~,idx_b] = min(D(idx_clus)) ;
        idx_A(a) = idx_clus(idx_b) ;
    end
else
    % delete the column one by one
    idx_A = ones(n,1) ;
    for a = n-1:-1:k
        % (a+1) columns, a clusters
        idx = kmeans(A(find(idx_A),:),a) ;
        count_a = zeros(a,1) ;
        % find the cluster including 2 columns
        for temp = 1:a+1
            count_a(idx(temp)) = count_a(idx(temp)) + 1 ;
        end
        clus_2 = find(count_a==2) ;
        % delete one colum in this cluster randomly
        idx_2 = find(idx == clus_2) ;
        idx_nonzero = find(idx_A) ;
        idx_2 = idx_nonzero(idx_2) ;
        if norm(A_save(:,idx_2(1)),2) < norm(A_save(:,idx_2(2)),2)
            idx_A(idx_2(1)) = 0 ;
        else
            idx_A(idx_2(2)) = 0 ;
        end
%         temp_rand = binornd(1,0.5) ;
%         idx_A(idx_2(temp_rand+1)) = 0 ;
%         A(idx_2(temp_rand+1),:) = [] ;
    end
    idx_A = find(idx_A) ;
end
% compute the projection matrix
A = A_save(:,idx_A) ;
Ps = A * pinv(A) ;

        
        
