alpha = 1.5;
sigma = 0.5;
c = 1.2;
tau = 1.5;
T = 1e-8;
n = 4000;
t = 0.1;

%% check psi
S = 0;
for i = 1:n
    W = Model2rnd(alpha, sigma, c, tau, T);
    S = S + exp(-t*sum(W));
end
fprintf(1, 'true %f, empirical %f\n', ...
    exp(-Model2psi(t, alpha, sigma, c, tau)), S/n);

%% check kappa
S = 0;
m = 3;
for i = 1:n
    W = Model2rnd(alpha, sigma, c, tau, T);
    S = S + sum(W.^m.*exp(-t*W));
end
fprintf(1, 'true %f, empirical %f\n', ...
    Model2kappa(m, t, alpha, sigma, c, tau), S/n);
    
