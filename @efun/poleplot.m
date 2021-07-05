function poleplot(s, varargin)
% poleplot(s) plots the poles of r(x) = ift(s) in the complex plane. 
% and superimposes the values of r(x) (plotted in the x-y plane) over the
% plot. 
%
% poleplot(s, 'zt') plots the subset of poles of modulus < 1 of the function r(z) = r(x), 
% where z = exp(2 * pi *1i *x). 
%
% See also: poles, plot. 

dom = s.domain; 
ms = 25;
res = 2*s.res+1; 
hold_state = ishold;

if ~isempty(varargin) && any(strcmp(varargin{:}, 'zt'))
    x = linspace(0, 2*pi, 1000); 
    plot(exp(1i*x), '-k', 'linewidth', 1.5)
    pol = poles(s, 'zt'); pol = pol(abs(pol)< 1); 
    w = s.weights; 
    plot(exp(linspace(0, 1,100)*2*pi*1i),'k', 'linewidth', 1.5)
    hold on
    plot(real(pol), imag(pol), '.', 'markersize', ms)
    axis equal
else %do the trig plot. 
    x = linspace(dom(1),dom(2),res);
    [H, ~, ~] = ift(s, 'polres'); 
    pol = poles(s, 'trig'); 
    plot(x, H(x), 'linewidth', 2, 'color','k')
    hold on
    plot(real(pol), imag(pol), '.b', 'markersize', ms)
    plot(dom, [0 0], '--k')  
end
if ~hold_state
    hold off
end

end


