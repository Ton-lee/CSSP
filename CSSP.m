% main program that conduct column subset selecting on a given matrix

m = 150;
n = 20;
k = 10;
rng(1);
A = randn(m,n);


% [idx_A, Ps] = kmeans_appro(A,k,2);  % k-means approach 
[idx_A, Ps] = traversal_appro(A,k);  % traversal approach

A_choose = A(:,idx_A) ;
A_approx = A_choose * pinv(A_choose) * A ;
[error, error_rate] = error_test(A, idx_A);
disp('A: ');
disp(A);
disp('Column choice: ');
disp(idx_A);
disp('Approximation: ');
disp(A_approx);
disp('Error (F-norm): ');
disp(error);
disp('Error (relative)');
disp(error_rate);
