function out = Model1augloglikel(cnts, ty, u, eta, sigma, c, tau)
n = sum(cnts);
y = c./(1 + exp(-ty));
out = (n-1)*log(u) - Model1psi(u, eta, sigma, c, tau) - gammaln(n);
out = out + length(cnts)*(log(eta) - gammaln(1-sigma) - log(c));
out = out + sum(gammaln(cnts-sigma) + ...
    (tau-sigma).*log(y) + log(c-y) - (cnts-sigma).*log(y + u));
end