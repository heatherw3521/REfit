function out = mtimes(r, a)
% a*m = times r by a, where a is a scalar. 
%
%%

if isempty(r)
    out = []; 
    return
elseif (isa(a, 'double') && isa(r, 'rfun')) || (isa(r, 'double') && isa(a, 'rfun'))
    out = times(r, a); 
else 
    error('efun:mtimes: mtimes is supported only with scalars.')
end


end