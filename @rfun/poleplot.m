function poleplot(r, varargin)
% plots the poles of r in the complex plane and superimposes the values 
% of r at a sample x (plotted in the x-y plane) over the plot. 
%
% poleplot(r, 'zt') plots the subset of poles of modulus < 1 of the 
% function r(z), where z = exp(2 * pi *1i *(x-r.domain(1))/diff(r.domain)). 
%
% See also: rfun/poles, rfun/plot. 

dom = r.domain; 
ms = 25; 
hold_state = ishold;

if ~isempty(varargin) && any(strcmp(varargin{:}, 'zt'))
    x = linspace(0, 2*pi, 1000); 
    plot(exp(1i*x), '-k', 'linewidth', 1.5)
    pol = poles(r, 'zt'); pol = pol(abs(pol)< 1);  
    hold on
    plot(real(pol), imag(pol), '.', 'markersize', ms)
    axis equal
else %plot on strip containing domain. 
    [s, x] = sample(r); 
    plot(x,s, 'linewidth', 2, 'color','k')
    hold on
    plot(real(r.poles), imag(r.poles), '.b', 'markersize', ms)
    plot(dom, [0 0], '--k')  
end

if ~hold_state
    hold off
end

end


