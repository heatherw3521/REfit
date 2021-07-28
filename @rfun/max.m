function mx = max(r, varargin)
% max(r) returns locations where local minima of r occur.
% max(r, 'abs') returns the location of the absolute maximum. 
%
%%
% see also rfun/extrema

[~,mx] = extrema(r); 

if ~isempty(varargin) %want absolute min
   [~, idx] = sort(feval(r, mx)); 
    mx = mx(idx(end)); 
end