clc;
close all;
clear all;
set(0,'defaultAxesFontSize', 16)
set(0,'DefaultTextFontSize', 20);
set(0, 'DefaultLineLineWidth', 1.5);
addpath(genpath('../utils'));
addpath(genpath('../GGP'));
addpath(genpath('..'));

%%
dataname = 'nips1000';
[cnts,words] = import_data(dataname);
cnts = cnts';
n = sum(cnts);

%%
n_samples = 50000;
results = GGPfit(cnts, n_samples);

%% posterior computation and prediction
results.pred.prop = [];
results.ks = [];
pi_post = [];
rank_temp = {};
thin = 100;
for i = 1:thin:length(results.eta)
    fprintf(1, '%d\n', i);
    W_pred = GGPrnd(results.eta(i), results.sigma(i), 1.0, 1e-6);
    cnts_pred = csizernd(W_pred, n);
    [~, ~, prop] = plot_loglog(cnts_pred, ...
        'xmin', min(cnts), 'xmax', max(cnts), 'show', false);
    results.pred.prop = [results.pred.prop; prop];
    rank_temp = [rank_temp; sort(cnts_pred, 'descend')];
    results.ks = [results.ks, compute_ks(cnts, cnts_pred)];
    W_post = gamrnd(cnts-results.sigma(i), 1./(results.u(i)+1));
    W_post_rem = GGPsumrnd(results.eta(i), results.sigma(i), results.u(i)+1);
    pi_post = [pi_post; W_post / (sum(W_post) + W_post_rem)];
end
maxrank = min(cellfun(@(x)length(x), rank_temp));
results.pred.rank = zeros(length(rank_temp), maxrank);
for i = 1 : length(rank_temp)
    results.pred.rank(i,:) = rank_temp{i}(1:maxrank);
end

results.post.pi_mean = mean(pi_post, 1);
results.post.pi_cred = quantile(pi_post, [0.025, 0.975], 1);

clear rank_temp;
clear pi_post;
clear W_post;
clear W_post_rem;
clear W_pred;
clear W_star;

%%
main_path = pwd();
chain_name = datestr(datetime('now'));
save_path = ['../results/' dataname '/GGP/chain ' chain_name];
mkdir(save_path);
cd(save_path);
save(['posterior ', chain_name]);




%%
figure;
plot(results.ll);
title('loglikelihood');
saveas(gcf,'loglikelihood.jpg');

figure;
histogram(results.eta, 50);
title('eta');
saveas(gcf,'posterior_eta.jpg');

figure;
histogram(results.sigma, 50);
title('sigma');
saveas(gcf,'posterior_sigma.jpg');

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
saveas(gcf,'pred_proportion.jpg');

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
saveas(gcf,'pred_rank.jpg');

%%
figure;
hold on
n_words = 10;
title('posteriors');
[sorted_cnts, ind] = sort(cnts, 'descend');
obs_W = sorted_cnts/sum(sorted_cnts);
obs_W = obs_W(1:n_words);
cred = results.post.pi_cred(:,ind);
for i = 1:n_words
    plot([i, i], [cred(1,i), cred(2,i)], 'r', 'linewidth', 3);
    plot(i, obs_W(i), 'g*');
end
xlabel('Rank')
ylabel('Prob')
xlim([1, n_words]);
legend({'95% credible interval posterior', 'Data'})
saveas(gcf,'posterior_pi.jpg');


%%
figure;
hold on
n_words = 50;
title('posteriors');
[sorted_cnts, ind] = sort(cnts, 'descend');
obs_W = sorted_cnts/sum(sorted_cnts);
obs_W = obs_W(1:n_words);
cred = results.post.pi_cred(:,ind);
for i = 1:n_words
    plot([i, i], [cred(1,i), cred(2,i)], 'r', 'linewidth', 3);
    plot(i, obs_W(i), 'g*');
end
xlabel('Rank')
ylabel('Prob')
xlim([1, n_words]);
legend({'95% credible interval posterior', 'Data'})
saveas(gcf,'posterior_pi_50.jpg');

cd(main_path);
