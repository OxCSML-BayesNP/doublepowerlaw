function W = Model2rnd(eta, sigma, c, tau, T)
eta0 = eta * exp(gammaln(tau) - tau*log(c));
W0 = GGPrnd(eta0, sigma, 1, T);
G = gamrnd(tau, 1./c, size(W0));
W = W0./G;
end