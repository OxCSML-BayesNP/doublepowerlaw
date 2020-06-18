function [ty, u, eta, sigma, c, tau] = Model2init(cnts)
eta = gamrnd(1000., 1);
sigma = sum(cnts==1)/length(cnts);
c = 1.0;
[~, tau] = estimate_changepoint(cnts);
S = sum(Model2rnd(eta, sigma, c, tau, 1e-8));
u = gamrnd(sum(cnts), 1/S);
ty = randn(size(cnts));
end