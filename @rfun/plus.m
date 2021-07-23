 function h = plus(r, g)
 % Compression plus: Add together r+g and represent the result as an rfun. 
 % One of r or g must be an rfun. The other can be a scalar, an efun, 
 % or an rfun. 
 %
 % plus(r, g) is called for synax 'r+g'. 
 %
 % See also efun/compress. 
  
%%
if ~isa(r, 'rfun') || isempty(r) % make the first input an rfun: 
    h = plus(g,r); 
    return

elseif isempty(g) %efun + []
    h = []; 
    return
    
elseif isa(g, 'double') %g is a scalar
    if ~all(size(g)==1)
        error('rfun:plus:input must be rfuns, scalars, or efuns.')
    end
    h = r; 
    h.const = r.const + g; % just add in the constant. 
    return
    
else
    % either g is an rfun or an efun
    if ~(r.domain == g.domain) %check that domains match
    error('rfun:plus: domain of definition for each object must be the same')
    end
    %check for simple case where g is double r. 
    if isa(g, 'rfun') && length(g) == length(r)
        if  all(abs( sort(abs(r.nodes))-sort(abs(g.nodes))) < 1e-15) &&...
                all(abs(sort(abs(r.weights)) - sort(abs(g.weights))) < 1e-15)
            h.const = 2*r.const; %double constant.
            h.vals = 2*r.vals; %double values at nodes. 
            return
        end
    end
    
    % call compression + in efun 
    % (this avoids using the rfun constructor directly and 
    % running into issues with sampling grids that don't account
    % for singularity locations):
    rf = ft(r); 
    if isa(g, 'rfun') 
      gf = ft(g); 
    else
      gf = g;
    end
    hf = rf+gf; 
    h = ift(hf); 
    %to do: add a 'type' flag allowing the return of an efun if desired. 
end
       
end