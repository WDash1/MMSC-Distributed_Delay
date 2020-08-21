%This function may be used to calculate the zeros and integral weights
%associated with the Trapezium rule approximation of an integral
%over a given interval.
function [zeros, weights] = computeTrapeziumRuleWeights(a, b, n)
    if(n<2)
        error('The n value must be >=2.');
    end
    
    zeros = linspace(a,b,n);
    h = (b-a)/(n-1);
    weights = ones(1,n);
    weights(1) = 0.5;
    weights(end) = 0.5;
    weights = h.*weights;
end