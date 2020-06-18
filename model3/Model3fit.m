function results = Model3fit(cnts, n_samples)

% initialize parameters
[u, eta, sigma, c, beta, w0, tau] = Model3init(cnts);
log_eta = log(eta);
logit_sigma = log(sigma) - log(1-sigma);
delta = tau - sigma;
log_delta = log(delta);
log_beta = log(beta);
log_w0 = log(w0);
log_u = log(u);

fprintf(1, 'initial estimates: eta %f, sigma %f, beta %f, w0 %f, tau %f\n',...
    eta, sigma, beta, w0, tau);

results.ll = [];
results.eta = [];
results.sigma = [];
results.beta = [];
results.w0 = [];
results.tau = [];
results.u = [];
burn_in = 0.5*n_samples;
thin = 100;

for i = 1:n_samples
    % resample eta
    log_eta_new = log_eta + 0.01*randn;
    eta_new = exp(log_eta_new);
    log_r = Model3loglikel(cnts, u, eta_new, sigma, c, beta, w0, tau) ...
        - 0.5*log_eta_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*log_eta^2;
    if rand < exp(min(log_r, 0.0))
        log_eta = log_eta_new;
        eta = eta_new;
    end

    % resample sigma
    logit_sigma_new = logit_sigma + 0.01*randn;
    sigma_new = 1/(1+exp(-logit_sigma_new));
    tau_new = sigma_new + delta;
    log_r = Model3loglikel(cnts, u, eta, sigma_new, c, beta, w0, tau_new) ...
        - 0.5*logit_sigma_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*logit_sigma^2;
    if rand < exp(min(log_r, 0.0))
        logit_sigma = logit_sigma_new;
        sigma = sigma_new;
        tau = tau_new;
    end

    % resample beta
    log_beta_new = log_beta + 0.01*randn;
    beta_new = exp(log_beta_new);
    log_r = Model3loglikel(cnts, u, eta, sigma, c, beta_new, w0, tau) ...
        - 0.5*log_beta_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*log_beta^2;
    if rand < exp(min(log_r, 0.0))
        log_beta = log_beta_new;
        beta = beta_new;
    end

    % resample delta
    log_delta_new = log_delta + 0.01*randn;
    delta_new = exp(log_delta_new);
    tau_new = sigma + delta_new;
    log_r = Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau_new) ...
        - 0.5*log_delta_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*log_delta^2;
    if rand < exp(min(log_r, 0.0))
        log_delta = log_delta_new;
        delta = delta_new;
        tau = tau_new;
    end
    
    % resample w0
    log_w0_new = log_w0 + 0.05*randn;
    w0_new = exp(log_w0_new);
    log_r = Model3loglikel(cnts, u, eta, sigma, c, beta, w0_new, tau) ...
        - 0.5*log_w0_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*log_w0^2;
    if rand < exp(min(log_r, 0.0))
        log_w0 = log_w0_new;
        w0 = w0_new;
    end
    
    % resample u
    log_u_new = log_u + 0.01*randn;
    u_new = exp(log_u_new);
    log_r = Model3loglikel(cnts, u_new, eta, sigma, c, beta, w0, tau) ...
        - 0.5*log_u_new^2 ...
        - Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau) ...
        + 0.5*log_u^2;
    if rand < exp(min(log_r, 0.0))
        log_u = log_u_new;
        u = u_new;
    end
    
    ll = Model3loglikel(cnts, u, eta, sigma, c, beta, w0, tau);   
    if i > burn_in
        results.ll = [results.ll, ll];
        results.eta = [results.eta, eta];
        results.sigma = [results.sigma, sigma];
        results.beta = [results.beta, beta];
        results.w0 = [results.w0, w0];
        results.tau = [results.tau, tau];
        results.u = [results.u, u];                
    end
    
    if mod(i,thin) == 0
        fprintf(1, 'iter %d, ll %.5e, eta %f, sigma %f, beta %f, w0 %f, tau %f, u %5e\n', ...
            i, ll, eta, sigma, beta, w0, tau, u);
    end        
end
end