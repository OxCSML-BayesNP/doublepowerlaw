function out = GGPaugloglikel(cnts, u, eta, sigma, c)
n = sum(cnts);
out = (n-1)*log(u) - GGPpsi(u, eta, sigma, c) - gammaln(n);
out = out + length(cnts)*(log(eta) - gammaln(1-sigma));
out = out + sum(gammaln(cnts-sigma) - (cnts-sigma)*log(u + c));
end
