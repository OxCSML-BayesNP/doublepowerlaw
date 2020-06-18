%%
clc;
close all;
clear all;
set(0,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',20);
set(0, 'DefaultLineLineWidth', 1.5);
addpath(genpath('../utils'));
addpath(genpath('../GGP'));

%%
eta_true = 1000.;
sigma_true = 0.2;
c = 1.0;
beta_true = 1000.;
w0_true = 1.0;
tau_true = 2.0;
T = 1e-8;
n = 500000;
[W_true, K_true] = Model3rnd(eta_true, sigma_true, c, beta_true, w0_true, tau_true, T);
cnts = csizernd(W_true, n);

%%
figure;
plot_loglog(cnts);
title('observed data - proportion');
figure;
plot_rank(cnts);
title('observed data - rank');

%%
n_samples = 7000;
results = Model3fit(cnts, n_samples);

%% prediction
results.pred.prop = [];
rank_temp = {};
thin = 100;
for i = 1:thin:length(results.eta)
    fprintf(1, '%d\n', i);
    W_pred = Model3rnd(results.eta(i), results.sigma(i), 1.0, ...
        results.beta(i), results.w0(i), results.tau(i), 1e-8);
    cnts_pred = csizernd(W_pred, n);
    [~, ~, prop] = plot_loglog(cnts_pred, ...
        'xmin', min(cnts), 'xmax', max(cnts), 'show', false);
    results.pred.prop = [results.pred.prop; prop];
    rank_temp = [rank_temp; sort(cnts_pred, 'descend')];
end
maxrank = min(cellfun(@(x)length(x), rank_temp));
results.pred.rank = zeros(length(rank_temp), maxrank);
for i = 1 : length(rank_temp)
    results.pred.rank(i,:) = rank_temp{i}(1:maxrank);
end
clear rank_temp;

%%
figure;
plot(results.ll);
title('loglikelihood');

figure;
histogram(results.eta, 50);
title('eta');

figure;
histogram(results.beta, 50);
title('beta');

figure;
histogram(results.sigma, 50);
title('sigma');

figure;
histogram(results.w0, 50);
title('w0');

figure;
histogram(results.tau, 50);
title('tau');
%%
figure;
title('prediction - proportion');
[~, centerbins, ~] = plot_loglog(cnts, 'show', false);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
hold on
plot_variance = @(x,lower,upper,color) fill([x,x(end:-1:1)], ...
    [upper,lower(end:-1:1)],color, 'EdgeColor', color);
quantile_prop = quantile(results.pred.prop, [.025, .975]);    
plot(centerbins, quantile_prop, 'color', [.8, .8, 1], 'linewidth', 2.5);
hold on
ind2 = quantile_prop(1,:)>0;
ha = plot_variance(centerbins(ind2), ...
    quantile_prop(1,ind2),quantile_prop(2,ind2), [.8, .8, 1] );
[~, centerbins, ~] = plot_loglog(cnts, 'linespec', 'ro');

%%
figure;
title('prediction - rank');
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
hold on
quantile_rank = quantile(results.pred.rank, [.025, .975]);
plot(1:maxrank, quantile_rank, 'color', [.8, .8, 1], 'linewidth', 2.5);
plot_variance(1:maxrank, quantile_rank(1,:), quantile_rank(2,:), [.8, .8, 1]);
plot_rank(cnts, 'color', 'r');

%%
save(['simulation ', datestr(datetime('now'))]);
