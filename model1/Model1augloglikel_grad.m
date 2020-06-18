function out = Model1augloglikel_grad(ty, cnts, u, sigma, c, tau)
y = c./(1 + exp(-ty));
out = ((tau - sigma)./y - 1./(c-y) - (cnts-sigma)./(y+u)).*y.*(c-y)./c;
end