function results = Model2fit(cnts, n_samples)

% initialize parameters
[ty, u, eta, sigma, ~, tau] = Model2init(cnts);
c = 1;
log_u = log(u);
log_eta = log(eta);
logit_sigma = log(sigma) - log(1-sigma);
delta = tau - sigma;
log_delta = log(delta);

fprintf(1, 'initial estimates: eta %f, sigma %f, tau %f\n',...
    eta, sigma, tau);

results.ll = [];
results.eta = [];
results.sigma = [];
results.tau = [];
results.u = [];

% Store variables for posterior pi
results.thinned_ll = [];
results.thinned_eta = [];
results.thinned_sigma = [];
results.thinned_tau = [];
results.thinned_u = [];
results.thinned_y = [];

burn_in = 0.5*n_samples;
thin = 100;
eps = 0.05;
L = 30;
for i = 1:n_samples
    % resample ty
    U = @(x)(-Model2augloglikel(cnts, x, u, eta, sigma, c, tau));
    grad_U = @(x)(-Model2augloglikel_grad(x, cnts, u, sigma, c, tau));
    ty = hmc_step(ty, U, grad_U, eps, L);
   
    % resample eta
    log_eta_new = log_eta + 0.01*randn;
    eta_new = exp(log_eta_new);
    log_r = Model2augloglikel(cnts, ty, u, eta_new, sigma, c, tau) ...
        - 0.5*log_eta_new^2 ...
        - Model2augloglikel(cnts, ty, u, eta, sigma, c, tau) ...
        + 0.5*log_eta^2;
    if rand < exp(min(log_r, 0.0))
        log_eta = log_eta_new;
        eta = eta_new;
    end

    % resample sigma
    logit_sigma_new = logit_sigma + 0.01*randn;
    sigma_new = 1/(1+exp(-logit_sigma_new));
    tau_new = sigma_new + delta;
    log_r = Model2augloglikel(cnts, ty, u, eta, sigma_new, c, tau_new) ...
        - 0.5*logit_sigma_new^2 ...
        - Model2augloglikel(cnts, ty, u, eta, sigma, c, tau) ...
        + 0.5*logit_sigma^2;
    if rand < exp(min(log_r, 0.0))
        logit_sigma = logit_sigma_new;
        sigma = sigma_new;
        tau = tau_new;
    end

    % resample delta
    log_delta_new = log_delta + 0.01*randn;
    delta_new = exp(log_delta_new);
    tau_new = sigma + delta_new;
    log_r = Model2augloglikel(cnts, ty, u, eta, sigma, c, tau_new) ...
        - 0.5*log_delta_new^2 ...
        - Model2augloglikel(cnts, ty, u, eta, sigma, c, tau) ...
        + 0.5*log_delta^2;
    if rand < exp(min(log_r, 0.0))
        log_delta = log_delta_new;
        delta = delta_new;
        tau = tau_new;
    end

    % resample u
    log_u_new = log_u + 0.01*randn;
    u_new = exp(log_u_new);
    log_r = Model2augloglikel(cnts, ty, u_new, eta, sigma, c, sigma+delta) ...
        - 0.5*log_u_new^2 ...
        - Model2augloglikel(cnts, ty, u, eta, sigma, c, sigma+delta) ...
        + 0.5*log_u^2;
    if rand < exp(min(log_r, 0.0))
        log_u = log_u_new;
        u = u_new;
    end
    
    ll = Model2augloglikel(cnts, ty, u, eta, sigma, c, tau);
    if i > burn_in
        results.ll = [results.ll, ll];
        results.eta = [results.eta, eta];
        results.sigma = [results.sigma, sigma];
        results.tau = [results.tau, tau];
        results.u = [results.u, u];                
    end
    
    if mod(i,thin) == 0
        fprintf(1, 'iter %d, ll %.5e, eta %f, sigma %f, tau %f, u %5e\n', ...
            i, ll, eta, sigma, tau, u);
    end        
    
    if (i > burn_in) && (mod(i,thin) == 0)
        results.thinned_ll = [results.thinned_ll, ll];
        results.thinned_eta = [results.thinned_eta, eta];
        results.thinned_sigma = [results.thinned_sigma, sigma];
        results.thinned_tau = [results.thinned_tau, tau];
        results.thinned_u = [results.thinned_u, u];
        results.thinned_y = [results.thinned_y; exp(ty)];
    end
end
end