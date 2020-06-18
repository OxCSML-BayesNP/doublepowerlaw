function [cp, tau] = estimate_changepoint(cnts)
figure;
sorted = sort(cnts, 'descend');
loglog(sorted);
hold on;
mdl = fitlm(log(100:200), log(sorted(100:200)));
beta = mdl.Coefficients.Estimate;
x0 = 0:.1:15;
y0 = beta(1) + beta(2)*x0;
loglog(exp(x0), exp(y0));

figure;
x = log(100:length(cnts));
y = beta(1) + beta(2)*x;
absdiff = abs(sorted(100:end)-exp(y));
loglog(exp(x), absdiff);
cp = 100 + sum(absdiff < 5);
tau = -1/beta(2);
end