function mins = min(s, varargin)
% find locations for all local minima of ift(s) on its domain.
% min(s, 'abs') returns a location for the absolute minimum. 
% 
% 
%%
% see also rfun/extrema,rfun/roots.

r = ift(s); 
mins = min(r, varargin{:}); 
end