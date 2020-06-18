function S = Model3sumrnd(eta, sigma, c, beta, w0, tau)
S = GGPsumrnd(eta, sigma, c);
K = poissrnd(beta*w0^(-tau)/tau);
S = S + sum(w0*rand(1, K).^(-1/tau));
end