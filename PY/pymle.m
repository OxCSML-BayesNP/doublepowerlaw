function [xmle, fval, exitflag] = pymle(m)

% find MLE of (theta,sigma) parameters for a PY partition

fun =@(x) -pypartitionpdf(m, x(1), x(2));
x0 = [1; sum(m==1)/length(m)]; % initialization
[xmle,fval,exitflag] = fminsearch(fun,x0);
fval= -fval;

end