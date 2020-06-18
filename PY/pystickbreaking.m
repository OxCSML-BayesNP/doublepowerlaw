function weights = pystickbreaking(alpha, sigma, K)

% Sample K first size-biased weights of a Pitman-Yor process using
% stick-breaking

b = betarnd(1-sigma, alpha + (1:K)*sigma);
logb = log(b);
log1mb = [0, cumsum(log(1.0-b(1:end-1)))];
logweights = logb + log1mb;
weights = exp(logweights);
