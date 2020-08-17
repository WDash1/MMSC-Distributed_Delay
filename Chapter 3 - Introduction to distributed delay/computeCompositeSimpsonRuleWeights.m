%This function may be used to calculate the zeros and integral weights
%associated with the composite Simpson rule approximation of an integral
%over a given interval.
function [zeros, weights] = computeCompositeSimpsonRuleWeights(a, b, n)
    zeros = linspace(a,b,n);
    
    if(n<3)
        error('The n value must be >=3');
    end
    
    if(mod(n,2) == 0)
        error('The n value must be odd.');
    end
    
    
    weights = ones(1,n).*2;
    weights((1:((n-mod(n,2))/2)).*2) = 4;
    weights(1) = 1;
    weights(end) = 1;
    
    h = (b-a)/(n-1);
    
    weights = (h/3) .* weights;
       
end