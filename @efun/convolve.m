function h = convolve(s,g)
% Compute the convolution of s and g (the product of the
% trigonometric rationals they represent.) 
% 
% s must be an efun. g can be an rfun or a function handle for a function 
% defined in Fourier space. If g is an rfun, then convolve(s,g) 
% computes convolve(s, ft(g)). 
%
% See also: rfun/times.

if ~isa(s, 'efun') % ?? times efun
    error('efun:convolve: first input argument must be an efun.')
    
elseif isempty(g) % efun times []
    h = efun([]); 
    return
       
else  %rfuns, efuns, function handles
    if isa(g, 'rfun') %if rfun, convert to efun
        g = ft(g); 
    end

%% do convolution as a product in value space: 
% check domains: 
    h = s; 
    if ~isa(g, 'function_handle')
        if ~all(s.domain == g.domain)
            error('efun:times: domain of definition for each must be the same.')
        end
        if any(strcmpi([g.space, s.space], 'time')) 
            error('efun:times: "convolve" for efuns in signal space not yet available.')
        end
        D = max(s.res, g.res); 
        m = min((D-1 +mod(D-1,2))/2, length(s)+length(g)); %max possible # terms in sum. 
        h.tol = max(s.tol, g.tol);
    else
        D = s.res;
        m = (D-1 +mod(D-1,2))/2;  
        h.tol = s.tol; 
    end
    %sample and take prod in time: 
    S = feval(s, (0:D).'); 
    G = feval(g, (0:D).'); 
    S = coeffs2sample(S, 'pos'); 
    G = coeffs2sample(G, 'pos'); 
    HH = sample2coeffs(S.*G, 'pos'); 
    %now use compression algorithm to construct h = conv(s,g):       
    tol = 1e-11; 
    happy = false;   
    M = min(2*m+20, D); 
    H = HH;
    h.const = HH(1); 
    H(1) = 0; 
    h.scl = max(abs(H)); 
    H = H/h.scl; 
    while ~happy
        coeffs = H(1:M+1); 
        [r,ss,w] = prony_compress(coeffs, m, tol); 
        h.exp = r; 
        h.weights = w;  
        err = abs(feval(h,(0:D).')- HH);
        if max(err) < tol && M<=D
            happy = true;
        elseif M == D
            break; 
        else
            M = min(2*M, D); %double sample # and try again
        end
    end
    h.sv = ss;
    h.res = D; 
end

end




