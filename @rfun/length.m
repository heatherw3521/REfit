function out = length( r )
% LENGTH OF AN RFUN. If r is a type (m-1, m) rational function, 
% then length(r) = m. 
%   
%%

% Empty check: 
if (isempty(r)) 
   out = []; 
   return
else
   out = length(r.nodes); 
end
