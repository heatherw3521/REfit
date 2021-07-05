function out = ldivide(s, g)
% s.\a = division of s by a, where a is a scalar. 
%
% 
%%
if isempty(s)
    out = []; 
elseif isa(g, double) && isa(s, efun)
    if g ==0 
        out = inf; 
    else 
    out = times(s, 1/a); 
    end
else 
    error('efun:ldivide: only division by scalars is supported.')
end


end