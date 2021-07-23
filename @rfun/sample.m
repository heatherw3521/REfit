function varargout = sample(r, varargin)
% sample an rfun. 
% [vals, loc] = sample(r) returns an equispaced sample (vals) of r on
% its domain. 
%
% See also: rfun/feval, rfun/getcoeffs
%%

if isempty(r)
    if nargout==2
        varargout = {[],[]}; 
    else
        varargout = {[]};
    end
    return
end

res = r.res; 
dom = r.domain; 
loc = linspace(dom(1), dom(2), res+1); loc = loc(1:end-1).';
vals = feval(r, loc); 

if nargout==2
    varargout = {vals, loc}; 
else
    varargout = {vals};
end

end
    
    