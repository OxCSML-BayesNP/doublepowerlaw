function W = Model1rnd(eta, sigma, c, tau, T)
eta0 = eta*c^(tau-sigma)/tau;
W0 = GGPrnd(eta0, sigma, c, T);
B = betarnd(tau, 1, size(W0));
W = W0./B;
end