%%
clc;
clear;

%%
a = [0.1; 0.3];
A = [2.25, 1.2;1.2, 1.25];
iA = inv(A);
U = @(x)(0.5*(x-a)'*iA*(x-a));
grad_U = @(x)(iA*(x-a));

%%
eps = 0.1;
L = 20;
N = 5000;
x = zeros(2, N);
x(:,1) = randn(2, 1);
for i = 2:N
    x(:,i) = hmc_step(x(:,i-1), U, grad_U, eps, L);
end

%%
error_ellipse('C', A, 'mu', a);
hold on;
scatter(x(1,:), x(2,:));