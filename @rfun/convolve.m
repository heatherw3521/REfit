function h = convolve(r,g, varargin)
% Compute the convolution between r and g.  
% One of r, g must be an rfun. The other can be a function handle, 
% or an rfun.
%
% One can also input r as an rfun and g as an efun. In this case, 
% h = convolve(r, ift(g)). 
%
% See also: efun/times. 

if ~isa(r, 'rfun') % ?? convolve rfun
    h = times(g,r); 
    return
    
elseif isempty(g) % rfun convolve []
    h = rfun([]); 
    return  
else
    rs = ft(r); 
    if isa(g, 'function_handle') % convolve r with a function
        [~, x] = sample(r); 
        gs = efun(g(x), 'domain', r.domain, 'tol', r.tol); %represent function with efun. 
    elseif isa(g, 'rfun') && all(g.domain == r.domain)
        gs = ft(g); 
    elseif isa(g, 'efun') && all(g.domain == r.domain)
        gs = g; 
    else
        error('rfun:convolve:domains must be the same')
    end  
hs = times(rs, gs);
h = ift(hs); 
end

