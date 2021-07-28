function mxs = max(s, varargin)
% find locations for all local mamxima of ift(s) on its domain.
% max(s, 'abs') returns a location for the absolute maximum. 
% 
% 
%%
% see also rfun/extrema,rfun/roots.

r = ift(s); 
mxs = max(r, varargin{:}); 
end