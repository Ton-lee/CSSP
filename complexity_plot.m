methods = {'kmeans1', 'kmeans2', 'greedy1', 'greedy2', 'deterministic', 'randomized'};
n = 2:2:64; m = 64; k=n/2; N = 20;

times = zeros(length(methods), length(n));
times(1,:) = (N+1)*k.*n;
times(2,:) = N/3*(n.^3-k.^3);
times(3,:) = 1/4*n.*k^4-1/30*k.^5;
times(4,:) = 1/10*(n.^5-k.^5);
times(5,:) = m*n.^2;
times(6,:) = 40*m*n.*min(m,n);