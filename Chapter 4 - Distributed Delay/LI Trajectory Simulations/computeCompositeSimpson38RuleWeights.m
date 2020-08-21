%This function may be used to calculate the zeros and integral weights
%associated with the composite Simpson 3/8 rule approximation of an integral
%over a given interval.
function [zeros, weights] = computeCompositeSimpson38RuleWeights(a, b, n)
    zeros = linspace(a,b,n);
    
    if(n<4)
        error('The n value must be >=4');
    end
    
    if(mod(n,3) ~= 1)
        error('The n value must be exactly one more than a multiple of 3.');
    end
    
    
    weights = ones(1,n).*3;
    weights(((1:((n-mod(n,3))/3)).*3)+1) = 2;
    weights(1) = 1;
    weights(end) = 1;
    
    h = (b-a)/(n-1);
    
    weights = (3*h/8) .* weights;
       
end