 function h = plus(s, g)
 % Compression plus: Add together s+g and represent the result as an efun. 
 % One of s or g must be an efun. The other can be a scalar, an efun, 
 % or an rfun. 
 %
 % plus(s, g) is called for synax 's+g'. 
 %
 %%
 % See also efun/compress. 
  
%%
if ~isa(s, 'efun') || isempty(s) % make the first input an efun: 
    h = plus(g,s); 
    return

elseif isempty(g) %efun + []
    h = []; 
    return
    
elseif isa(g, 'double') %g is a scalar
    if ~all(size(g)==1)
        error('efun:plus:input must be efuns, scalars, or rfuns.')
    end
    h = s; 
    h.const = s.const+g; %just add the constant value in. 
    return
    
else  %rfuns and efuns
    if isa(g, 'rfun') %convert to efun
        g = ft(g); 
    end
    if any(strcmpi([g.space, s.space], 'time') )%still need to set up 
                                               % compression for signal
                                               % space.
        error('efun:plus:"+" for efuns in signal space not yet available.')
    end
    h = s; 
    if ~(s.domain == g.domain) %check that domains match
    error('efun:times: domain of definition for each must be the same')
    end
    %check for simple case where g is double s. 
    if length(g) == length(s) && g.const == s.const
        if  all(abs(s.exp-g.exp) < 1e-15) &&...
                all(abs(s.scl*s.weights - g.scl*g.weights) < 1e-15)
            h.scl = 2*s.scl; %just double
            h.const = 2*s.const; 
            return
        end
    end
    %  deal with constants and scales. 
    s.weights = s.scl*s.weights; %push scl into weights. 
    g.weights = g.scl*g.weights; 
    const = s.const + g.const; 
    h.scl = 1; s.scl = 1; g.scl = 1; %reset scls and consts. 
    h.const = 0; s.const = 0; g.const = 0; 

    % add two efuns with compression algorithm        
    tol = 1e-10; 
    happy = false;
    D = max([s.res, g.res]); 
    MM = min( (length(s.exp)+ length(g.exp)), (D + mod(D,2))/2);
    M = min(2*MM+20, D); 
    coeffs = feval(s,(0:M).') + feval(g,(0:M).'); 
    scl = max(abs(coeffs)); 
    if scl==0
        %efun is interpolating zero function
        h.exp = -.5; %this is arbitrary
        h.weights = 0; 
        h.const = 0; 
        h.scl = 1; 
        h.res = 1; 
        h.tol = tol;
        return
    end
        
    coeffs = coeffs/scl; 
    while ~happy
        [r,ss,w] = prony_compress(coeffs, MM, tol); 
        h.exp = r; 
        h.weights = w;
        h.scl = scl; 
        rnd = randi([0 D], 1,50).'; 
        err = abs(feval(h,[(0:M).';rnd])-[coeffs; feval(s,rnd) + feval(g,rnd)]); 
        if max(err) < tol 
            happy = true;
        end
        if M > 2*D
            break; 
        end
        coeffs = [scl*coeffs; feval(s, (M+1:2*M).')+feval(g, (M+1:2*M).')];
        scl = max(abs(coeffs)); 
        coeffs = coeffs/scl; 
        M = 2*M; %double sample # and try again
    end
    h.sv = ss;
    h.scl = scl; 
    h.const = const; 
    h.res = D; 
end

end