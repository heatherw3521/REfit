function out = integral(s, varargin)
% integration of an efun over a fixed interval. 
%
% integral(s) integrates ift(s) over its domain. 
%
% integral(s, [a, b]) integrates ift(s) over the interval [a, b]. 
%
% See also, sum, cumsum. 

out = sum(s, varargin{:}); 
end