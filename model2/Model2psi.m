function y = Model2psi(t, eta, sigma, c, tau)
y = (eta/sigma) * integral(...
    @(x)(((x+t).^sigma - x.^sigma).* x.^(tau-sigma-1) .* exp(-c*x)), 0, inf);
end