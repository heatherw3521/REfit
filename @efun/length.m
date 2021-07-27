function out = length( F )
%length(F) = m.  F is a sum m of weighted exponentials.
%   
%%

% Empty check: 
if (isempty(F)) 
   out = []; 
   return
else
   out = length(F.exp); 
end
