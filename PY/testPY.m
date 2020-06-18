close all
clear all


%% Sample from PY
theta = 100;
sigma = 0.5;
n = 100000;
[~, m, K, ~] = pycrprnd(theta, sigma, n);

figure
loglog(sort(m, 'descend'))
xlabel('Rank')
ylabel('Counts')

%% Compute partition for range of values of sigma
sigma_all = [max(0, -theta):.001:.9];
for i=1:length(sigma_all)
    [logproba1(i)] = pypartitionpdf(m, theta, sigma_all(i));
end

figure
plot(sigma_all, exp(logproba1-max(logproba1)))
hold on; plot(sigma, 0, '*g')
xlabel('sigma')
ylabel('likelihood')

%% Compute partition for a range of values of theta
theta_all = -sigma+.001:.1:200;
for i=1:length(theta_all)
    [logproba2(i)] = pypartitionpdf(m, theta_all(i), sigma);
end

figure
plot(theta_all, exp(logproba2 - max(logproba2)))
hold on; plot(theta, 0, '*g')
xlabel('theta')
ylabel('likelihood')

%% Find MLE by optimization
[xmle,fval,exitflag] = pymle(m);

theta_mle = xmle(1)
sigma_mle = xmle(2)
% Check it has found the right value
pypartitionpdf(m, theta, sigma)
pypartitionpdf(m, xmle(1), xmle(2))


%% ANC data

%% Load data and find MLE of PY parameters
load ../data/nips/nips1000.mat;
wordcounts = cnts;
[xmle,fval,exitflag] = pymle(wordcounts);


% Check it makes sense
theta_mle = xmle(1)
sigma_mle = xmle(2)
sum(wordcounts ==1)/length(wordcounts)

% Look at the predictive for large counts
n = sum(wordcounts);
nsamples = 1000;
for i=1:nsamples
    % Use approximation: otherwise n is too large to sample exactly
    out(i, :) = pylargecountsapproxrnd(theta_mle, sigma_mle, n, 5);
end
wordcounts_sort = sort(wordcounts, 'descend');

figure
for i=1:5
    plot([i, i], quantile(out(:, i), [.025, .975]), 'r', 'linewidth', 3)
    hold on
    plot(i, wordcounts_sort(i), '*g')
end
xlim([0.5, 5.5])
xlabel('Rank')
ylabel('Counts')
legend({'95% credible interval predictive', 'Data'})

% Sanity check: check that  if data come from exact PY model with fitted parameters, all is
% fine
ntest = 100000;
[~, wordcounts2, K, ~] = pycrprnd(theta_mle, sigma_mle, ntest);
nsamples = 1000;
for i=1:nsamples
    out(i, :) = pylargecountsapproxrnd(theta_mle, sigma_mle, ntest, 5);
end

wordcounts2_sort = sort(wordcounts2, 'descend');
figure
for i=1:5
    plot([i, i], quantile(out(:, i), [.025, .975]), 'r', 'linewidth', 3)
    hold on
    plot(i, wordcounts2_sort(i), '*g')
end
xlim([0.5, 5.5])
xlabel('Rank')
ylabel('Counts')
legend({'95% credible interval predictive', 'Data'})
