function out = sum(r, varargin)
% integration of an rfun over a fixed interval. 
% sum(r) integrates r over its domain. 
% sum(r, [a, b]) integrates r over the interval [a, b]. 
%%
% See also, rfun/cumsum, rfun/integral. 

%%
if isempty(r)
    out = []; 
    return
end

if nargin == 1
    %integral over domain
       rdom = r.domain; 
       out = (rdom(2)-rdom(1))*r.const; 
else %integral over a fixed interval
    dom = varargin{1}; 
    a = dom(1); 
    b = dom(2); 
    rdom = r.domain; 
    if a < rdom(1) || rdom(2) < b
        error('rfun:sum:bounds of integral must be a subset of [%d %d]',rdom(1), rdom(2)); 
    end
    h = cumsum(r);
    out = h(b) -h(a); 
end

end


