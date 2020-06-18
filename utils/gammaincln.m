function y = gammaincln(x, a)
y = zeros(size(a));
neg = a < 0;
if sum(neg > 0)
    if length(x) > 1
        xn = x(neg);
    else
        xn = x;
    end
    an = a(neg);
    y(neg) = log(exp(gammaincln(xn, an+1)) - xn.^an.*exp(-xn)) - log(an);
end
zero = a == 0;
if sum(zero)
    if length(x) > 1
        xz = x(zero);
    else
        xz = x;
    end
    y(zero) = log(expint(xz));
end
pos = a > 0;
if sum(pos)
    if length(x) > 1
        xp = x(pos);
    else
        xp = x;
    end
    ap = a(pos);
    y(pos) = gammaln(ap) + log(gammainc(xp, ap, 'upper'));
end
y = real(y);
