function [reLambda] = DispersRel(tau, a, b)

[f,g] = funs(tau,a,b,30);

values = roots(f,g,'resultant');
real_componnents = values(:,1);

reLambda = max(real_componnents);


end

