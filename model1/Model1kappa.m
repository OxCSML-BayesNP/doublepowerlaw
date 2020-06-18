function y = Model1kappa(m, t, eta, sigma, c, tau)
y = eta * exp(gammaln(m-sigma) - gammaln(1-sigma));
y = y * integral(@(x)((x+t).^(sigma-m).*x.^(tau-sigma-1)), 0, c);
end