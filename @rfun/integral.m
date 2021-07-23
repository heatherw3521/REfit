function out = integral(r, varargin)
% integration of an rfun over a fixed interval. 
%
% integral(r) integrates r over its domain. 
% integral(r, [a, b]) integrates ift(r) over the interval [a, b]. 
%
% See also, rfun/sum, rfun/cumsum. 

out = sum(r, varargin{:}); 
end