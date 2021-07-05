function out = length( F )
%LENGTH OF AN EFUN = m, where F represents a length m sum of weighted exponentials.
%   
%%

% Empty check: 
if (isempty(F)) 
   out = []; 
   return
else
   out = length(F.exp); 
end
