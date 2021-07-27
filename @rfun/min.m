function mins = min(r, varargin)
% min(r) returns locations where local minima of r occur.
% min(r, 'abs') returns the location of the absolute minimum. 
%
%%
% see also rfun/extrema

[mins, ~, ~] = extrema(r); 

if ~isempty(varargin) %want absolute min
   [~, idx] = sort(feval(r, mins));
   mins = mins(idx(1)); 
end
