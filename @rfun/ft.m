function s = ft(r, varargin)
% Compute the Fourier transform of an rfun. 
% This returns an efun object (sum of complex exponentials)
% representing the Fourier transform of the trigonometric rational r.
%
% For details, see: 
%
% See also: efun/efun, efun/ift

%%
% check for tol parameter: 
tol = 1e-10;
if ~isempty(varargin)
    tolid = find(strcmpi(varargin, 'tol'));
    tol = varargin{tolid+1};
end

poles = r.poles;  
dom = r.domain; 
a = dom(1); b = dom(2); 
%%
% the exponents are known in closed form via the poles of r.  
% We average the conjugate pairs of poles since conjugate symmetry is not 
% enforced by rfun:

p = exp(-2i*pi/(b-a)*poles); 
p1 = p(abs(p)<1); 
p2 = p(abs(p)>1); 
p1 = sort(p1, 'ComparisonMethod', 'real');
p2 = 1./conj(p2); 
p2 = sort(p2, 'ComparisonMethod', 'real');

% average the poles: 
p = (p1 + p2)/2;
%poles = log(p)/(-2i*pi/(b-a)); 
%p = exp(-2i*pi/(b-a)*poles); 
s = efun(); 
s.nodes = log(p);
%%
% get a sample in Fourier space to solve for weights:
res = r.res;
res = res + (1-mod(res, 2));
happy = false; 
while ~happy
    x = linspace(a, b, res+2); 
    x = x(1:end-1);
    y = get_coeffs_exp(feval(r,x), 'pos'); %does an fft
    %check for resolution: (check decay on coeffs):
        %y = y((resolve+1):end);
    cutoff = standardChop(y, tol);
    happy = cutoff < length(y); 
    res = 2*res; 
end
yy = y(1:cutoff); 
%%
% use a subset of the samples to solve for weights:

mm = length(poles); 
happy = false; 
while ~happy 
    mm = min(2*mm, cutoff); 
    %mm = length(yy); 
    y = yy(1:mm); 
    %%
    % solve for the weights: 
    x = (0:mm-1).';
    n = length(p); 
    V = repmat((p.'), mm,1); 
    V = V.^(repmat( x, 1,n));  
    w = V\y; 
    %check for small weights: 
    if ~isempty(find(abs(w)<= tol, 1))
        p = p((abs(w)> tol)); 
        s.nodes = log(p);  
        n = length(p); 
        %need to recompute the weights 
        V = repmat((p.'), mm,1); 
        V = V.^(repmat( x, 1,n));  
        w = V\y;
    end
    %check accuracy: 
    if mm < cutoff
        sp = min(cutoff, 20); 
        xx = randi([0, cutoff-1], sp,1); 
        nn = length(xx); 
        V = repmat((p.'), nn,1); 
        V = V.^(repmat(xx, 1,length(p)));  
        happy = norm(V*w-yy(xx+1)) < 1e1*tol;
    else
        happy =true; 
    end
end

%%
% construct efun
scl = norm(w);
w = w/scl; 
s.weights = w; 
s.domain = r.domain; 
s.space = 'Fourier'; 
s.domain = [-(cutoff-1) cutoff+1]; 
s.const = r.const; 
s.scl = scl;

end






