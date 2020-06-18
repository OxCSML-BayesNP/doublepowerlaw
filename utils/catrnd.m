function x = catrnd(p, n)
if nargin < 2
    n = 1;
end
u = rand(1, n) * sum(p);
bins = [0, cumsum(p)];
[~, x] = histc(u, bins);
end