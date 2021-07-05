function plotcoeffs(s, varargin)
%
% semilogy plot of the magnitude of the Fourier coefficients 
% associated with r(x) = ift(s). 
%
%%

if isempty(s)
    return %nothing to do. 
end

res = s.res; %resolution parameter. 
coeffs = abs(feval(s, 0:res-1)); 
holdstate = ishold; 
ms = 25; 
%% build plot
semilogy((-(res-1):(res-1)).', [flip(coeffs(2:end)); coeffs], '.', 'markersize', ms)
% set some properties: 
if ~holdstate
    set(gcf,'color','w')
    set(gca, 'fontsize', 18)
    title('Fourier coefficients')
    hold off    
end

end




