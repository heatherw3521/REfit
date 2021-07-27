function s = uminus(s)
%  unary minus.
%  -s negates the efun s.
%
%   uminus(s) is called for the syntax '-s'.
%%

if ( isempty(s) )
    return
end

s.scl = -s.scl; 
s.const = -s.const; 
end