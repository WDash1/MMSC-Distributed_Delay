function [f,g] = funs(tau,a,b,dom)
u = a+b;v = b/(a+b)^2;
%Chebfun2 area to look for roots in
d = [-dom, dom, -0.1, dom];
%Write f and g as the real and imaginary parts of the dispersion relation,
%where lambda=x+iy for real x,y. The full dispersion relation is given by:
%lambda^2 + lambda+lambda*(u^2-2*u*v)*exp(-lambda*tau)+u^2*exp(-lambda*tau)


%characteristic = @(z)   ((u.^2 + z) .* (z+1+4.*u.*v - 6.*u.*v.*exp(-z.*tau))) ...
%                        -2.*u.*v.*(2.*u - 3.*u.^2.*exp(-z.*tau));

%characteristic = @(z)   z.^2 + z.*(1+4.*u.*v+u.^2 - 6.*u.*v.*exp(-z.*tau)) ...
%                        +u.^2 - 4.*u.^2.*v + 4.*u.^3.*v;
                    


f = chebfun2(@(x,y)x.^2-y.^2 + x.*(1+4.*u.*v + u.^2 - 6.*u.*v.*exp(-x.*tau).*cos(y.*tau))...,
                    -y.*(6.*u.*v.*exp(-x.*tau).*sin(y.*tau))...,
                    +u.^2 - 4.*u.^2.*v+4.*u.^3.*v, d);
g = chebfun2(@(x,y)2.*x.*y + y.*(1+4.*u.*v + u.^2 - 6.*u.*v.*exp(-x.*tau).*cos(y.*tau))...,
                    +x.*(6.*u.*v.*exp(-x.*tau).*sin(y.*tau)),d);
%v((u.^2 + z) .* (z+1+4.*u.*v - 6.*u.*v.*exp(-z.*tau))) ...
%                        -2.*u.*v.*(2.*u - 3.*u.^2.*exp(-z.*tau));



%characteristic = @(z) z.^2 + z.*(1+(u.^2 - 2.*u.*v).*exp(-z.*tau)) + u.^2.*exp(-z.*tau);


%f = chebfun2(@(x,y)real(characteristic(x+y.*1j)), d);
%g = chebfun2(@(x,y)imag(characteristic(x+y.*1j)), d);

%f = chebfun2(@(x,y)x.^2-y.^2+x...,
%+(x.*cos(tau*y)+y.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%+u^2*exp(-x*tau).*cos(tau*y),d);
%g = chebfun2(@(x,y)2*x.*y +y...,
%+(y.*cos(tau*y)-x.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%-u^2*exp(-x*tau).*sin(tau*y),d);
end

