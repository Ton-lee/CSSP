% svd approximation approach
function A_appro = svd_appro(A,k)
[U,D,V] = svd(A) ;
A_appro = U(:,1:k) * D(1:k,1:k) * V(:,1:k)' ;
