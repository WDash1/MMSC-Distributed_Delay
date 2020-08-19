function [f,g] = funs(tau,a,b,dom)
u = a+b;v = b/(a+b)^2;
%Chebfun2 area to look for roots in
d = [-dom, dom, -0.1, dom];
%Write f and g as the real and imaginary parts of the dispersion relation,
%where lambda=x+iy for real x,y. The full dispersion relation is given by:
%lambda^2 + lambda+lambda*(u^2-2*u*v)*exp(-lambda*tau)+u^2*exp(-lambda*tau)

characteristic = @(z) z.^2 + z.*(1+(u.^2 - 2.*u.*v).*exp(-z.*tau)) + u.^2.*exp(-z.*tau);


f = chebfun2(@(x,y)real(characteristic(x+y.*1j)), d);
g = chebfun2(@(x,y)imag(characteristic(x+y.*1j)), d);


%f = chebfun2(@(x,y)x.^2-y.^2+x...,
%+(x.*cos(tau*y)+y.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%+u^2*exp(-x*tau).*cos(tau*y),d);
%g = chebfun2(@(x,y)2*x.*y +y...,
%+(y.*cos(tau*y)-x.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%-u^2*exp(-x*tau).*sin(tau*y),d);
end

