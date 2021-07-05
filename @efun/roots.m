function z = roots(s, varargin)
% Find the roots of r(x) = ift(s), the inverse Fourier transform of s. 
% 
% roots(s) finds the real-valued roots of r(x) on its domain. 
%
% roots(s, 'all') finds all roots, including complex valued ones. 
%
% See also: rfun/roots

%%
% call rfun to compute roots: 
if isempty(s)
    z = []; 
else
    z = roots(ift(s), varargin{:}); 
end
