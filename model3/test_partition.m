%%
close all;
clear;
clc;

%% sample CRM from model 3
eta = 1000.0;
sigma = 0.2;
c = 1.0;
beta = 1000.0;
w0 = 1.0;
tau = 2.00;
T = 1e-10;
[W, K] = Model3rnd(eta, sigma, c, beta, w0, tau, T);
fprintf(1, 'W has %d atoms (%d core atoms)\n', length(W), K);

%% sample partitions
n_list = [10^4, 10^5, 10^6, 10^7, 10^8];
for i = 1:length(n_list)
    csize{i} = csizernd(W, n_list(i));
    leg{i} = ['n = ' num2str(n_list(i))];
end

%% plot proportions of cluster sizes
figure
for i = 1:length(n_list)
    plot_loglog(csize{i});
    hold on;
end
legend(leg)
hold off;

%% plot ranks
figure
xlabel('Rank')
ylabel('Frequency')
for i = 1:length(n_list)
    plot_rank(csize{i});
    hold on;
end
legend(leg);
hold off;