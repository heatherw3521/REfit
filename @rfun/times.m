function h = times(r,g)
% Compute the product of r and g.  
% One of r, g must be an rfun. The other can be a scalar or an rfun.
%
% times(r, g) is called when the syntax 'r.*g' is used. 
%
% See also: efun/convolve. 


%%
tol = 1e-10; 

if ~isa(r, 'rfun') % ?? times rfun
    h = times(g,r); 
    return
    
elseif isempty(g) % rfun times []
    h = rfun([]); 
    return
    
elseif isa(g, 'double') %g is a scalar
    if ~all(size(g)==1)
        error('efun:times:input must be rfuns, scalars, or efuns.')
    end
    h = r; 
    h.vals = g*r.vals; 
    h.const = g*r.const; 
    return
    
else  %rfuns and efuns
    %check domains: 
    if ~all(s.domain == g.domain)
       error('efun:times: domain of definition for each must be the same')
    end
    % to compute r.*g,  we convert to Fourier space and convolve
    % the sums. This helps with the issue of accounting for 
    % singularity locations. (Auto-sampling/simple sampling schemes
    % in the rfun constructor can be an issue here; efun seems to do better.)  
    
    rf = ft(r); 
    if isa(g, 'rfun') %convert to efun
        gf = ft(g); 
    end
    hf = times(rf, gf); 
    h = ift(hf); 
end
