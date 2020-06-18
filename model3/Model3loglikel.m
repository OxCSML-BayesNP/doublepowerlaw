function y = Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau)
n = sum(cnts);
y = (n-1)*log(u) -Model3psi(u, eta, sigma, c, beta, w0, tau) - gammaln(n) ...
    + sum(Model3logkappa(cnts, u, eta, sigma, c, beta, w0, tau));
end