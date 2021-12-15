function [idx_A,A_appro] = greedy_appro(A,k,mode) 
% mode 1 : greedily add
% mode 2 : greedily cut
% output : the index and the approximation of A
[~,n] = size(A) ;
if mode == 1
    % initialize an empty choice
    idx_choice = zeros(n,1) ;
    A_choose = [] ;
    for a = 1:k
        % find the candiate choice
        idx_candiate = find(idx_choice == 0) ;
        temp_error = zeros(length(idx_candiate),1) ;
        for b = 1:length(idx_candiate)
            A_temp = [A_choose,A(:,idx_candiate(b))] ;
            A_appro = A_temp * pinv(A_temp) * A ;
            temp_error(b) = norm(A-A_appro,'fro') ;
        end
        % find the best one to add
        [~,idx_now] = min(temp_error) ;
        idx_choice(idx_candiate(idx_now)) = 1 ;
        % update the choice
        A_choose = [A_choose,A(:,idx_now)] ;
    end
    idx_A = find(idx_choice==1) ;
    Ps = A_choose * pinv(A_choose) ;
elseif mode == 2
    % initialize a full choice
    idx_choice = ones(n,1) ;
    for a = 1:n-k
        idx_candiate = find(idx_choice==1) ;
        temp_error = zeros(length(idx_candiate),1) ;
        for b = 1:length(idx_candiate)
            idx_temp = [idx_candiate(1:b-1);idx_candiate(b+1:end)];
            A_temp = A(:,idx_temp) ;
            A_appro = A_temp * pinv(A_temp) * A ;
            temp_error(b) = norm(A-A_appro,'fro') ;
        end
        % find the one with limited influncence to cut
        [~,idx_now] = min(temp_error) ;
        idx_choice(idx_candiate(idx_now)) = 0 ;
    end
    idx_A = find(idx_choice == 1) ;
    A_choose = A(:,idx_A) ;
    Ps = A_choose * pinv(A_choose) ;
end
A_appro = Ps * A ;

    
                
                