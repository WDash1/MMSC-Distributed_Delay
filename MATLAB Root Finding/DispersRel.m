function [reLambda] = DispersRel(tau, a, b)

[f,g] = funs(tau,a,b,20);

lambda = max(roots([f;g]));
if(abs(lambda(1))<1e-10)
    lambda(1)=-1;
end

reLambda = lambda(1);
%imLambda = lambda(2);

end

