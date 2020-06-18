function x = hmc_step(x0, U, grad_U, eps, L)
x = x0;
r = randn(size(x));
H0 = U(x) + 0.5*sum(r.*r);
for i = 1:L
    r = r - 0.5*eps*grad_U(x);
    x = x + eps*r;
    r = r - 0.5*eps*grad_U(x);
end
H = U(x) + 0.5*sum(r.*r);
if rand > exp(min(H0-H, 0))
    x = x0;
end
end