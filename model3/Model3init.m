function [u, eta, sigma, c, beta, w0, tau] = Model3init(cnts)
sigma = sum(cnts==1)/length(cnts);
c = 1.0;
w0 = 1.0;
[cp, tau] = estimate_changepoint(cnts);
beta = tau*cp*w0^tau;
w0 = 10.0;
eta = beta;
S = Model3sumrnd(eta, sigma, c, beta, w0, tau);
u = gamrnd(sum(cnts), 1/S);
end