function h = times(s,g)
% Compute the product of s and g. If s and g are
% exponential sums in Fourier space, this is the convolution of the
% trigonometric rationals they represent. 
% 
% One of s, g must be an efun. The other can be a scalar or an rfun.
%
% times(s, g) is called when the syntax 's.*g' is used. 
%
% See also: rfun\convolve.

tol = 1e-10; 

if ~isa(s, 'efun') % ?? times efun
    h = times(g,s); 
    return
    
elseif isempty(g) % efun times []
    h = efun([]); 
    return
    
elseif isa(g, 'double') %g is a scalar
    if ~all(size(g)==1)
        error('efun:times:input must be efuns, scalars, or rfuns.')
    end
    h = s; 
    h.scl = g*s.scl; 
    h.const = g*s.const; 
    return
    
else  %rfuns and efuns
    if isa(g, 'rfun') %convert to efun
        g = ft(g); 
    end
    if any(strcmpi([g.space, s.space], 'time')) %still need to set up 
                                               % compression for signal
                                               % space.
        error('efun:times: "times" for efuns in signal space not yet available.')
    end
%% do compression times: 

% part 1: guess the upper bound on length of h: 
m1 = length(s).*length(g); 

%check for insignificant weight combinations:
G = (s.weights)*g.weights.';  
G = G(:); 
[~, idxw] = find(abs(G) < tol); % # weight combos below tol

%check for small exponents (usually a more efficient way): 
G = s.exp + (g.exp).';
G = G(:); 
[~, idxe] = find(abs(G) > abs(log(tol))); %exponents too big.
out = union(idxw, idxe); 

%set the upper bound: 
m = m1 - length(out); 

%use compression algorithm to construct h = s.*g:       
    tol = 1e-10; 
    happy = false; 
    D = max([s.res, g.res]); 
    M = min(2*m+20, D); 
    coeffs = feval(s,0:M).*feval(g,0:M); 
    while ~happy
        h.const = coeffs(1); 
        h.scl = max(abs(coeffs)); 
        coeffs = (coeffs-coeffs(1))/scl;
        [r,ss,w] = prony_compress(coeffs, MM, tol); 
        h.exp = r; 
        h.weights = w;
        rnd = randi([0 D], 1,50).'; 
        err = abs(feval(h,[(0:M).';rnd])-([h.scl*coeffs+h.const;feval(s,rnd)]).*...
        [h.scl*coeffs+h.const;feval(g,rnd)]); 
        if max(err) < tol 
            happy = true;
        end
        if M > 2*D
            break; 
        end
        coeffs = [coeffs; feval(s, M+1:2*M).*feval(g, M+1:2*M)];
        M = 2*M; %double sample # and try again
    end
    h.svs = ss;
    h.res = D; 
end

end