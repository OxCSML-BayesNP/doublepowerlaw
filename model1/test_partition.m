%%
close all;
clear;
clc;

%% sample CRM from model 1
alpha = 1000.0;
sigma = 0.2;
c = 2.0;
tau = 3.0;
T = 1e-10;
W = Model1rnd(alpha, sigma, c, tau, T);
fprintf(1, 'W has %d atoms\n', length(W));

%% sample partitions
n_list = [10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10];
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