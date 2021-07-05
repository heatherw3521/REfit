function out = sum(s, varargin)
% integration of an efun over a fixed interval. 
%
% sum(s) integrates ift(s) over its domain. 
%
% sum(s, [a, b]) integrates ift(s) over the interval [a, b]. 
%
% See also, cumsum, integral. 

%%
% TO DO: allow the integration of the exponential sum over a
% finite interval in Fourier space. 

if isempty(s)
    out = []; 
    return
end

if length(varargin) == 1
    %integral over domain
       out = s.const; 
else %integral over a fixed interval
    dom = varargin{1}; 
    sdom = s.domain; 
    if dom(1) < sdom(1) || dom(2) > sdom(2)
        error('efun:sum:bounds of integral must be a subset of [%d %d]',sdom(1), sdom(2)); 
    end
    S = cumsum(s, dom(1), 'efun'); 
    out = S(dom(2)); 
end

end


