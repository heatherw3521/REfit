function h = times(r,g)
% Compute the product of r and g.  
% One of r, g must be an rfun. The other can be a scalar, rfun, 
% or function handle.
%
% One can also input r as an rfun and g as an efun. 
% Then, h =times(r,ift(g)).
%
% times(r, g) is called when the syntax 'r.*g' is used. 
%
% See also: efun/convolve. 

%%

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
    
elseif isa(g, 'function_handle')
    [rs, x] = sample(r); 
    h = rfun(rs.*g(x), x); 
    return
    
else  %rfuns and efuns
    %check domains: 
    if ~all(s.domain == g.domain)
       error('efun:times: domain of definition for each must be the same')
    end
    %for now we just call constructor. Maybe a better choice is to 
    % use convolve in efun instead:
    resr = r.res; 
    if isa(g, 'efun')
        resg = 2*(g.res)+1; 
    else
        resg = g.res; 
    end
    resm = max(resg, resr); 
    [rs, x] = sample(r, resm); 
    if isa(g, 'efun')
        gs = feval(g, 'values'); 
    else
        gs = feval(g); 
    h = rfun(gs.*rs, x); 
    
    %todo: add in the option to return efun. 
    %% alternative method:
    % to compute r.*g,  we convert to Fourier space and convolve
    % the sums. This can help with the issue of accounting for 
    % singularity locations. (Auto-sampling/simple sampling schemes
    % in the rfun constructor can be an issue here; efun seems to do better.)  
    
    %rf = ft(r); 
    %if isa(g, 'rfun') %convert to efun
    %    gf = ft(g); 
    %end
    %hf = convolve(rf, gf); 
    %h = ift(hf); 
end
