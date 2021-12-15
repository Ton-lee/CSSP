function H_random = lowrank_random(m,n,r)
% a random m *n matrix whose rank is r
A = randn(m,m) ;
B = randn(n,n) ;
[U,~,~] = svd(A) ;
[~,~,V] = svd(B) ;
s = sort(randn(r,1),'descend') ;
H_random = zeros(m,n) ;
for a = 1:r
    H_random = H_random + s(a) *power(1.2,-a) * U(:,a) * V(:,a).' ;
end
