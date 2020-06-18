function [u, eta, sigma, c] = GGPinit(cnts)
eta = gamrnd(1000., 1);
sigma = sum(cnts==1)/length(cnts);
c = 1.0;
S = GGPsumrnd(eta, sigma, c);
u = gamrnd(sum(cnts), 1/S);
end
