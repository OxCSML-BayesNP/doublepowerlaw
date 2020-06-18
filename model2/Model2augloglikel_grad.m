function out = Model2augloglikel_grad(ty, cnts, u, sigma, c, tau)
y = exp(ty);
out = tau - sigma - c.*y - (cnts-sigma).*y./(y+u);
end