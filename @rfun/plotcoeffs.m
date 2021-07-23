function plotcoeffs(r, varargin)
%
% semilogy plot of the magnitude of the Fourier coefficients 
% associated with r. 
%
%%

if isempty(r)
    return %nothing to do. 
end

s = ft(r); 
plotcoeffs(s)

end




