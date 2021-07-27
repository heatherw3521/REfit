function r = uminus(r)
%  unary minus.
%  -r negates the rfun r.
%
%   uminus(r) is called for the syntax '-r'.
%%

if ( isempty(r) )
    return
end

r.scl = -r.scl; 
r.const = -r.const; 
end