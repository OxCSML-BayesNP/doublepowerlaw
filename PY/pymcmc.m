function results = pymcmc(cnts, n_samples)
eta = gamrnd(1000., 1.);
sigma = sum(cnts==1)/sum(cnts);

log_eta = log(eta);
logit_sigma = log(sigma) - log(1-sigma);
fprintf(1, 'initial estimates: eta %f sigma %f\n', ...
    eta, sigma);

results.ll = [];
results.eta = [];
results.sigma = [];
burn_in = 0.5*n_samples;
thin = 100;
for i = 1:n_samples
    % resample eta
    log_eta_new = log_eta + 0.01*randn;
    eta_new = exp(log_eta_new);
    log_r = pypartitionpdf(cnts, eta_new, sigma) ...
        - 0.5*log_eta_new^2 ...
        - pypartitionpdf(cnts, eta, sigma) ...
        + 0.5*log_eta^2;
    if rand < exp(min(log_r, 0.0))
        log_eta = log_eta_new;
        eta = eta_new;
    end

    % resample sigma
    logit_sigma_new = logit_sigma + 0.01*randn;
    sigma_new = 1/(1+exp(-logit_sigma_new));
    log_r = pypartitionpdf(cnts, eta, sigma_new) ...
        - 0.5*logit_sigma_new^2 ...
        - pypartitionpdf(cnts, eta, sigma) ...
        + 0.5*logit_sigma^2;
    if rand < exp(min(log_r, 0.0))
        logit_sigma = logit_sigma_new;
        sigma = sigma_new;
    end

    ll = pypartitionpdf(cnts, eta, sigma);
    if i > burn_in
        results.ll = [results.ll, ll];
        results.eta = [results.eta, eta];
        results.sigma = [results.sigma, sigma];
    end

    if mod(i,thin) == 0
        fprintf(1, 'iter %d, ll %.5e, eta %f, sigma %f\n', ...
            i, ll, eta, sigma);
    end
end
end
