function out = sum(s, varargin)
% integration of an efun over a fixed interval. 
% sum(s) integrates ift(s) over its domain. 
% sum(s, [a, b]) integrates ift(s) over the interval [a, b]. 
%%
% See also, cumsum, integral. 

%%
if isempty(s)
    out = []; 
    return
end

if nargin == 1
    %integral over domain
       sdom = s.domain; 
       out = (sdom(2)-sdom(1))*s.const; 
else %integral over a fixed interval
    dom = varargin{1}; 
    a = dom(1); 
    b = dom(2); 
    sdom = s.domain; 
    if a < sdom(1) || sdom(2) < b
        error('efun:sum:bounds of integral must be a subset of [%d %d]',sdom(1), sdom(2)); 
    end
    S = cumsum(s, 'efun');
    vv = feval(S, dom.', 'values');
    out = vv(2) -vv(1) +(b-a)*s.const; 
end

end


