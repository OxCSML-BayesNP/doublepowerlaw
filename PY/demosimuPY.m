close all
clear all


%% Sample from PY

sigma = 0.5;
n = 100000;


theta_all =  [-0.4, -.3, -.2, -.1, 0, 1];
out = zeros(length(theta_all), n);
for i=1:length(theta_all)
    theta = theta_all(i);
    [~, m, K, ~] = pycrprnd(theta, sigma, n);
    out(i, 1:K) = m;
end


figure
loglog(sort(out, 2, 'descend')', 'o')
xlabel('Rank')
ylabel('Counts')
legend(num2str(theta_all))