function results = GGPfit(cnts, n_samples)

[u, eta, sigma, ~] = GGPinit(cnts);
% fix c
c = 1;
log_u = log(u);
log_eta = log(eta);
logit_sigma = log(sigma) - log(1-sigma);

fprintf(1, 'initial estimates: u %f eta %f sigma %f\n', ...
    u, eta, sigma);

results.ll = [];
results.eta = [];
results.sigma = [];
results.u = [];
burn_in = 0.5*n_samples;
thin = 100;
for i = 1:n_samples
    % resample eta
    log_eta_new = log_eta + 0.01*randn;
    eta_new = exp(log_eta_new);
    log_r = GGPaugloglikel(cnts, u, eta_new, sigma, c) ...
        - 0.5*log_eta_new^2 ...
        - GGPaugloglikel(cnts, u, eta, sigma, c) ...
        + 0.5*log_eta^2;
    if rand < exp(min(log_r, 0.0))
        log_eta = log_eta_new;
        eta = eta_new;
    end

    % resample sigma
    logit_sigma_new = logit_sigma + 0.01*randn;
    sigma_new = 1/(1+exp(-logit_sigma_new));
    log_r = GGPaugloglikel(cnts, u, eta, sigma_new, c) ...
        - 0.5*logit_sigma_new^2 ...
        - GGPaugloglikel(cnts, u, eta, sigma, c) ...
        + 0.5*logit_sigma^2;
    if rand < exp(min(log_r, 0.0))
        logit_sigma = logit_sigma_new;
        sigma = sigma_new;
    end

    % resample u
    log_u_new = log_u + 0.01*randn;
    u_new = exp(log_u_new);
    log_r = GGPaugloglikel(cnts, u_new, eta, sigma, c) ...
        - 0.5*log_u_new^2 ...
        - GGPaugloglikel(cnts, u, eta, sigma, c) ...
        + 0.5*log_u^2;
    if rand < exp(min(log_r, 0.0))
        log_u = log_u_new;
        u = u_new;
    end

    ll = GGPaugloglikel(cnts, u, eta, sigma, c);
    if i > burn_in
        results.ll = [results.ll, ll];
        results.eta = [results.eta, eta];
        results.sigma = [results.sigma, sigma];
        results.u = [results.u, u];
    end

    if mod(i,thin) == 0
        fprintf(1, 'iter %d, ll %.5e, eta %f, sigma %f, u %5e\n', ...
            i, ll, eta, sigma, u);
    end
end
end
