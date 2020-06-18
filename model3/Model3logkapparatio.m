function y = Model3logkapparatio(m, z, eta, sigma, c, beta, w0, tau)
loga = log(m-sigma) - log(z+c);
logb = log(beta) + (m-sigma)*log(z+c) - log(eta) - (m+1-tau)*log(z) ...
    + gammaln(1-sigma) + gammaincln(w0*z, m+1-tau) - gammaln(m-sigma);
logc = log(beta) + (m-sigma)*log(z+c) - log(eta) - (m-tau)*log(z) ...
    + gammaln(1-sigma) + gammaincln(w0*z, m-tau) - gammaln(m-sigma);
if length(loga) == 1
    y = logaddexp(loga,logb) - logaddexp(0, logc);
else
    if size(loga,1) == 1
        y = logsumexp([loga;logb],1) - logsumexp([zeros(size(logc));logc],1);
    else
        y = logsumexp([loga,logb],2) - logsumexp([zeros(size(logc)),logc],2);
    end
end