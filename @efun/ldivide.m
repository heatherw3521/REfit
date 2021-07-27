function out = ldivide(s, a)
% s.\a = division of s by a, where a is a scalar. 
%%
if isempty(s)
    out = []; 
elseif isa(a, double) && isa(s, efun)
    if a ==0 
        out = inf; 
    else 
    out = times(s, 1/a); 
    end
else 
    error('efun:ldivide: only division by scalars is supported.')
end


end