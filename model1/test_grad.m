%%
clc;
clear;

%%
alpha = 2.3;
sigma = 0.2;
c = 1.2;
tau = 3.2;
W = Model1rnd(alpha, sigma, c, tau, 1e-8);
n = 10000;
cnts = csizernd(W, n);

%%
u = 1000;
ty = randn(size(cnts));
grad = Model1augloglikel_grad(ty, cnts, u, sigma, c, tau);

%%
ll = Model1augloglikel(cnts, ty, u, alpha, sigma, c, tau);
ngrad = zeros(size(cnts));
for i = 1:length(cnts)
    dty = zeros(size(ty));
    dty(i) = 1e-5;
    llp = Model1augloglikel(cnts, ty+dty, u, alpha, sigma, c, tau);
    ngrad(i) = (llp - ll)/1e-5;
end

%%
fprintf(1, '%f\n', mean(abs(grad-ngrad)/abs(grad)));

    