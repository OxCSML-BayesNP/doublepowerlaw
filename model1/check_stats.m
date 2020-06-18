alpha = 1.5;
sigma = 0.5;
c = 1.2;
tau = 1.5;
T = 1e-8;
n = 5000;
t = 0.1;
%% check psi
S = 0;
for i = 1:n
    W = Model1rnd(alpha, sigma, c, tau, T);
    S = S + exp(-t*sum(W));
end
fprintf(1, 'true %f, empirical %f\n', ...
    exp(-Model1psi(t, alpha, sigma, c, tau)), S/n);

%% check kappa
S = 0;
m = 2;
for i = 1:n
    W = Model1rnd(alpha, sigma, c, tau, T);
    S = S + sum(W.^m.*exp(-t*W));
end
fprintf(1, 'true %f, empirical %f\n', ...
    Model1kappa(m, t, alpha, sigma, c, tau), S/n);
    
