function m = Model3urnrnd(n, eta, sigma, c, beta, w0, tau)
S = Model3sumrnd(eta, sigma, c, beta, w0, tau);
u = gamrnd(n, 1./S);
m = 1;
for i = 2:n
    K = length(m);
    if mod(i, 100) == 0
        fprintf(1, 'iter %d, %d clusters\n', i, K);
    end    
    p = Model3logkapparatio(m, u, eta, sigma, c, beta, w0, tau);
    p = [p, Model3logkappa(1, u, eta, sigma, c, beta, w0, tau)];
    p = exp(p - logsumexp(p, 2));
    k = catrnd(p);
    if k <= K        
        m(k) = m(k) + 1;
    else
        m(k) = 1;
    end
end
end