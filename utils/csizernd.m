function cnts = csizernd(W, n)
cnts = poissrnd(n*W/sum(W));
cnts = cnts(cnts>0);
end
