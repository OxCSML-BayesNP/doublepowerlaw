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
dataname='twitter_friends_300000';
n_chains = 3;

%% load chais
root = sprintf('../results/%s/PY', dataname);
d = dir(root);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
ws = cell(1, n_chains);
for i = 1:n_chains
    path = sprintf('%s/%s', root, nameFolds{i});
    Files = dir(fullfile(path,'*.mat'));    
    ws{i} = load(sprintf('%s/%s', path, Files.name));
end

%% collect all the samples   
cnts = ws{1}.cnts;
n = sum(cnts);
results.eta = [ws{1}.results.eta, ws{2}.results.eta, ws{3}.results.eta];
results.sigma = [ws{1}.results.sigma, ws{2}.results.sigma, ws{3}.results.sigma];
clear ws;

%% prediction
results.pred.prop = [];
results.ks = [];
rank_temp = {};
thin = 300;
T = 1e-6;
for i = 1:thin:length(results.eta)
    fprintf(1, '%d\n', i);
    W_pred = pystickbreaking(results.eta(i), results.sigma(i), 20*length(cnts));
    cnts_pred = csizernd(W_pred, n);
    [~, ~, prop] = plot_loglog(cnts_pred, ...
        'xmin', min(cnts), 'xmax', max(cnts), 'show', false);
    results.pred.prop = [results.pred.prop; prop];
    rank_temp = [rank_temp; sort(cnts_pred, 'descend')];
    results.ks = [results.ks, compute_ks(cnts, cnts_pred)];
end
maxrank = min(cellfun(@(x)length(x), rank_temp));
results.pred.rank = zeros(length(rank_temp), maxrank);
for i = 1 : length(rank_temp)
    results.pred.rank(i,:) = rank_temp{i}(1:maxrank);
end
clear rank_temp;

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
fprintf(1, 'KS: %f+-%f\n', mean(results.ks), std(results.ks));

%%
save(sprintf('%s/summary.mat', root), 'results');