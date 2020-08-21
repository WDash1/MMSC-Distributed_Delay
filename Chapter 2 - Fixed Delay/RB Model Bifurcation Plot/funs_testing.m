function [f,g] = funs_testing(tau,a,b,dom)
u = a+b;v = b/(a+b)^2;
%Chebfun2 area to look for roots in
d = [-dom, dom, -dom, dom];
%Write f and g as the real and imaginary parts of the dispersion relation,
%where lambda=x+iy for real x,y. The full dispersion relation is given by:
%lambda^2 + lambda+lambda*(u^2-2*u*v)*exp(-lambda*tau)+u^2*exp(-lambda*tau)


%characteristic = @(z)   ((u.^2 + z) .* (z+1+4.*u.*v - 6.*u.*v.*exp(-z.*tau))) ...
%                        -2.*u.*v.*(2.*u - 3.*u.^2.*exp(-z.*tau));




characteristic = @(z) z.^2 + z.*(1+(u.^2 - 2.*u.*v).*exp(-z.*tau)) + u.^2.*exp(-z.*tau);


f = chebfun2(@(x,y)real(characteristic(x+y.*1j)), d);
g = chebfun2(@(x,y)imag(characteristic(x+y.*1j)), d);

%sol = roots([f;g])









%f = chebfun2(@(x,y)x.^2-y.^2+x...,
%+(x.*cos(tau*y)+y.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%+u^2*exp(-x*tau).*cos(tau*y),d);
%g = chebfun2(@(x,y)2*x.*y +y...,
%+(y.*cos(tau*y)-x.*sin(tau*y)).*exp(-x*tau).*(u^2-2*u*v)...,
%-u^2*exp(-x*tau).*sin(tau*y),d);

sol = roots(f,g,'marchingsquares')
complex_values = sol(:,1) + sol(:,2).*1j
evaluation = characteristic(complex_values)



figure('Renderer', 'painters', 'Position', [10 10 300 300], 'Visible', 'on')
hold on
box on
plot(roots(f), 'b')
plot(roots(g), 'r')
plot(sol(:,1),sol(:,2),'g*')
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

set(gca,'XTick',d(1):50:d(2));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.0f'))
    
set(gca,'YTick',d(3):50:d(4));
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.0f'))

xlim([d(1) d(2)])
ylim([d(3) d(4)])

set(gca, 'FontSize', 15)
filename = './RB_Model_Roots_tau='+string(tau)+'_a='+string(a)+'_b='+string(b)+'.eps';
print('-depsc', '-tiff', '-r300', '-painters', filename);

end
