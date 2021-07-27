function [mins, maxs, infls] = extrema(r)
% find locations for all local minima/maxima/inflection points of r on its domain. 
% 
% 
%%
% see also rfun/min, rfun/max, rfun/inflection, rfun/roots.

% To find extrema, we differentiate, build an rfun, find its roots:
rd = diff(r,1, 'rfun'); 
ext = roots(rd); 
h = diff(rd); %handle for second derivative

%determine types of the extrema 
type = h(ext); 
mins = ext(type > 1e-8);
maxs = ext(type < 1e-8); 
infls = ext ( abs(type) < 1e-8); 

end
    



