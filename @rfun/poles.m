function p = poles(r, varargin)
% poles(r) returns the poles of r. 
%
% poles(r, 'zt') returns the poles of the rational r(z), 
% where z = exp(2*pi*1i*(x-r.domain(1))/diff(r.domain)). 
%
% See also: rfun/roots. 

if isempty(r)
    p = [];
    return
end

p = r.poles; 
if isempty(varargin)
    return
else
    dom = r.domain; 
    L = dom(2) - dom(1); 
    p = exp(2*pi*1i*(p-dom(1))/L); 
        %p = p(abs(p)<1); %only return poles inside disk.
end
    
end

