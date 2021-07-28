function [mins, maxs] = extrema(s)
% find locations for all local minima/maxima of ift(s) on its domain. 
% 
% 
%%
% see also rfun/extrema,rfun/roots.

r = ift(s); 
[mins, maxs] = extrema(r); 
end