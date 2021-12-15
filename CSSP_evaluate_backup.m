% test time and error performance under different k (with n fixed)
clear all; warning off;
% methods = {'blind', 'kmeans1', 'kmeans2', 'greedy1', 'greedy2', 'randomized', 'QR', 'svd', 'traversal'};
methods = {'blind', 'kmeans1', 'kmeans2', 'greedy1', 'greedy2', 'randomized', 'QR', 'svd'};
% methods = {'blind', 'kmeans1', 'kmeans2', 'QR', 'randomized'};
total_method = length(methods);
m = 64;
n = 32;
r = 24;
repeat = 100;
% k_list = 1:min(n-1, m);
k_list = 2:2:30;
time_list = zeros(total_method,length(k_list));
error_list = zeros(total_method, length(k_list));

for k_index = 1:length(k_list)
    k = k_list(k_index);
    disp('k='+string(k));
    time = zeros(total_method, repeat);
    error = zeros(total_method, repeat);
    for rep = 1:repeat
%         A = randn(m, n);
        A = randn(m, r) * randn(r, n);
        for method_i = 1:total_method
            tic;
            if strcmp(methods(method_i), 'svd')
                A_approx = svd_appro(A,k);
                time(method_i, rep) = toc;
                error_ = norm(A - A_approx, 'fro');
                error_rate = error_ / norm(A,'fro');
                error(method_i, rep) = error_rate;
            else
                if strcmp(methods(method_i), 'kmeans1')
                    [idx_A, Ps] = kmeans_appro(A,k,1);  % k-means approach
                elseif strcmp(methods(method_i), 'kmeans2')
                    [idx_A, Ps] = kmeans_appro(A,k,2);  % k-means approach
                elseif strcmp(methods(method_i), 'traversal')
                    [idx_A, Ps] = traversal_appro(A,k);  % traversal approach
                elseif strcmp(methods(method_i), 'blind')
                    [idx_A, Ps] = blind_appro(A,k);  % random approach
                elseif strcmp(methods(method_i), 'QR')
                    [idx_A, Ps] = QR_appro(A,k);  % QR approach
                elseif strcmp(methods(method_i), 'randomized')
                    [idx_A, Ps] = randomized_appro(A,k);  % randomized approach
                elseif strcmp(methods(method_i), 'greedy1')
                    [idx_A, Ps] = greedy_appro(A,k,1);  % greedy approach
                elseif strcmp(methods(method_i), 'greedy2')
                    [idx_A, Ps] = greedy_appro(A,k,2);  % greedy approach
    %             elseif strcmp(methods(method_i), 'svd')
    %                 [idx_A, Ps] = svd_appro(A,k);  % svd approach
                else
                    disp('Wrong method! (What is '+string(methods(method_i))+' ?)')
                end
                time(method_i, rep) = toc;
                [~, error_rate] = error_test(A, idx_A);
                error(method_i, rep) = error_rate;
            end
        end
    end
    time_list(:, k_index) = mean(time, 2);
    error_list(:, k_index) = mean(error, 2);
%     disp('.')
end
save(strcat('results/', datestr(now, 30), '.mat'));
%%
% figure;
% subplot(2,1,1)
% semilogy(k_list, time_list); title('time - k'); legend(methods);
% subplot(2,1,2);
% plot(k_list,10*log(error_list)); title('error - k'); legend(methods);
% plot(k_list,error_list); title('error - k'); legend(methods);

names = {
    'blind choose', ...
    'clustering center-based', ...
    'kmeans-based dropout', ...
    'greedy choice-based', ...
    'greedy dropout-based', ...
    'randomized', ...
    'deterministic', ...
    'svd', ...
    'traversal', ...
    };
method_names = names(1:length(methods));
figure;
plot(k_list, error_list(1,:), ':', ...
    k_list, error_list(2,:), '-', ...
    k_list, error_list(3,:), '-', ...
    k_list, error_list(4,:), '-', ...
    k_list, error_list(5,:), '-', ...
    k_list, error_list(6,:), '-', ...
    k_list, error_list(7,:), '-', ...
    k_list, error_list(8,:), '--' ...
    ); 
hold on;
if length(k_list) == 9
    plot( k_list, error_list(9,:), '--');
end
grid on; legend(method_names); xticks(k_list);