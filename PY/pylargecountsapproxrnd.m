function m = pylargecountsapproxrnd(theta, sigma, n, K)

% Sample the K's first counts for a Pitman-Yor process with parameter theta
% and sigma with n observations
% Uses a CLT approximation + stick-breaking

if n<10000 %
    warning('Large-scale approximation may be invalid')
end

weights = pystickbreaking(theta, sigma, K*100);% generates first weights using stick-breaking
m = randn(size(weights)).*sqrt(n*weights.*(1-weights)) + n*weights; % Use CLT approximation
m = sort(m, 'descend');
m = m(1:K);