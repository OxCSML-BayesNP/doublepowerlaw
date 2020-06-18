function weights = pypoststickbreaking(m, alpha, sigma)
% Sample K first size-biased weights of a posterior 
% Pitman-Yor process using stick-breaking

K = length(m);
b = betarnd(1-sigma + m, alpha + (1:K)*sigma + sum(m)-cumsum(m));
logb = log(b);
log1mb = [0, cumsum(log(1.0-b(1:end-1)))];
logweights = logb + log1mb;
weights = exp(logweights);

