function y = Model3logkappa(m, z, eta, sigma, c, beta, w0, tau)
loga = log(eta) + gammaln(m-sigma) - (m-sigma)*log(z+c) - gammaln(1-sigma);
logb = log(beta) + gammaincln(w0*z, m-tau) - (m-tau)*log(z);
y = logaddexp(loga, logb);
end