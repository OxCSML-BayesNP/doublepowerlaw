function [W, K] = Model3rnd(eta, sigma, c, beta, w0, tau, T)
[W, ~] = GGPrnd(eta, sigma, c, T);
K = poissrnd(beta*w0^(-tau)/tau);
Wcore = w0./(rand(1, K).^(1./tau));
W = [W, Wcore];
end