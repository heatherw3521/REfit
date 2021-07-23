function out = sum(r, varargin)
% integration of an rfun over a fixed interval. 
%
% sum(r) integrates r over its domain. 
% sum(r, [a, b]) integrates r over the interval [a, b]. 
%
% See also, cumsum, integral. 

%%
if isempty(r)
    out = []; 
    return
end

if length(varargin) == 1
    %integral over domain: since r is periodic, this will just be
    % the mean value of r over the domain * length of domain:
       dom = r.domain; 
       out = (dom(2)-dom(1))*r.const; 
else %integral over a fixed interval
    dom = varargin{1}; 
    rdom = r.domain; 
    if dom(1) < rdom(1) || dom(2) > rdom(2)
        error('efun:sum:bounds of integral must be a subset of [%d %d]',sdom(1), sdom(2)); 
    end
    S = cumsum(r, dom(1), 'efun'); %make a function handle
    out = S(dom(2)); %evaluate
end

end


