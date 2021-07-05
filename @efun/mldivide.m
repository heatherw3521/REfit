function out = mldivide(s, a)
% s\a = division of s by a scalar a. 
%
%%
if isempty(s)
    out = []; 
elseif a==0
    out = inf; 
else
    out = times(s, 1/a); 
end
    