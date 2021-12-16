% error testing function. We calculate the F-norm error as well as the relative error
function [error,error_rate] = error_test(A,index_num) 
% to compute the error caused by the approximation
% input : the initial matrix A & the choice of column index index_num
A_choose = A(:,index_num) ;
A_approx = A_choose * pinv(A_choose) * A ;
error = norm(A - A_approx, 'fro') ;
error_rate = error / norm(A,'fro') ;

