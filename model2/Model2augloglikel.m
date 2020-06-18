function out = Model2augloglikel(cnts, ty, u, eta, sigma, c, tau)
n = sum(cnts);
y = exp(ty);
out = (n-1)*log(u) - Model2psi(u, eta, sigma, c, tau) - gammaln(n);
out = out + length(cnts)*(log(eta) - gammaln(1-sigma));
out = out + sum(gammaln(cnts-sigma) + ...
    (tau-sigma).*log(y) - c.*y - (cnts-sigma).*log(y + u));
end