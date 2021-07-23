function out = ldivide(r, a)
% r.\a = division of r by a, where a is a scalar. 
%
% 
%%
if isempty(r)
    out = []; 
elseif isa(a, double) && isa(r, rfun)
    if a ==0 
        out = inf; 
    else 
    out = times(r, 1/a); 
    end
else 
    error('efun:ldivide: only division by scalars is supported.')
end


end