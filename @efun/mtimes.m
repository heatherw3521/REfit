function out = mtimes(s, a)
% a*m = times s by a, where a is a scalar. 
%
%%

if isempty(s)
    out = []; 
    return
elseif (isa(a, 'double') && isa(s, 'efun')) || (isa(s, 'double') && isa(a, 'efun'))
    out = times(s, a); 
else 
    error('efun:mtimes: mtimes is supported only with scalars.')
end


end