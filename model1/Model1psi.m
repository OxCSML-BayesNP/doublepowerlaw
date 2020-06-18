function y = Model1psi(t, eta, sigma, c, tau)
y = (eta/sigma) * integral(...
    @(x)(((t+x).^sigma - x.^sigma).*x.^(tau-sigma-1)), 0, c);
end