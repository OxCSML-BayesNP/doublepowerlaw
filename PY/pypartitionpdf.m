function [logproba, proba] = pypartitionpdf(m, theta, sigma)

% Computes the probability of the partition for a Pitman-Yor partition

% Check parameters
if (sigma<0) || (sigma>=1) || (theta <=-sigma)
    logproba = -Inf;
    proba = 0;
    return;
end

m = m(m>0);
K = length(m);
n = sum(m);
if sigma>0
    logproba = sum(gammaln(m-sigma)) - K*gammaln(1-sigma)...
        + gammaln(K + theta/sigma) + (K-1)* log(sigma) - gammaln(theta/sigma+1)...
        + gammaln(theta+1) - gammaln(theta + n);
else
    logproba = sum(gammaln(m))...
        + (K-1) * log(theta)...
        + gammaln(theta+1) - gammaln(theta + n);
end
proba = exp(logproba);
end