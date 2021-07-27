function plotcoeffs(r, varargin)
%
% semilogy plot of the magnitude of the Fourier coefficients 
% associated with r. 
%
%%

if isempty(r)
    return 
end

s = ft(r); 
plotcoeffs(s)

end




