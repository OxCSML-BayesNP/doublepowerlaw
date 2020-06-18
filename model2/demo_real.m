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
dataname = 'ANC';

[cnts,words] = import_data(dataname);
cnts = cnts';
n = sum(cnts);
c_true = 1;

%%
figure;
plot_loglog(cnts);
title('observed data - proportion');
figure;
plot_rank(cnts);
title('observed data - rank');

%%
n_samples = 100000;
results = Model2fit(cnts, n_samples);


%% prediction
results.pred.prop = [];
results.ks = [];
rank_temp = {};
thin = 100;
T = 1e-6;
for i = 1:thin:length(results.eta)
    W_pred = Model2rnd(results.eta(i), results.sigma(i), ...
        1.0, results.tau(i), T);
    cnts_pred = csizernd(W_pred, n);
    [~, ~, prop] = plot_loglog(cnts_pred, ...
        'xmin', min(cnts), 'xmax', max(cnts), 'show', false);
    results.pred.prop = [results.pred.prop; prop];
    results.ks = [results.ks, compute_ks(cnts, cnts_pred)];
    rank_temp = [rank_temp; sort(cnts_pred, 'descend')];
    W_post = gamrnd(cnts-results.sigma(i), 1./(results.u(i)+1));
end
maxrank = min(cellfun(@(x)length(x), rank_temp));
results.pred.rank = zeros(length(rank_temp), maxrank);
for i = 1 : length(rank_temp)
    results.pred.rank(i,:) = rank_temp{i}(1:maxrank);
end
clear rank_temp;


% posterior
pi_post = [];
for i=1:10:length(results.thinned_eta)
    curr_eta = results.thinned_eta(i);
    curr_sigma = results.thinned_sigma(i);
    curr_tau = results.thinned_tau(i);
    curr_u = results.thinned_u(i);
    curr_y = results.thinned_y(i,:);
    W_post = gamrnd(cnts-curr_sigma, 1./(curr_u+curr_y));
    W_star = Model2rnd(curr_eta, curr_sigma, c_true, curr_tau, T);
    z_ = binornd(1,exp(-curr_u*W_star));
    W_post_rem = sum(W_star.*z_);
    pi_post = [pi_post; W_post / (sum(W_post) + W_post_rem)];
end    

results.post.pi_mean = mean(pi_post, 1);
results.post.pi_cred = quantile(pi_post, [0.025, 0.975], 1);

clear pi_post;
fields = {'thinned_ll',...
            'thinned_eta',...
            'thinned_sigma',...
            'thinned_tau',...
            'thinned_u',...
            'thinned_y'};
results = rmfield(results,fields);
clear W_post;
clear W_post_rem;
clear W_pred;
clear W_star;
clear z_;

%%
main_path = pwd();
chain_name = datestr(datetime('now'));
save_path = ['../results/' dataname '/model2/' 'chain ' chain_name];
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

figure;
histogram(results.tau, 50);
title('tau');
saveas(gcf,'posterior_tau.jpg');

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
% Posterior
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


%%
% Posterior
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


cd(main_path);



