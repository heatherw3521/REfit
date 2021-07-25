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
res = []; 
if ~isempty(varargin)
    if nargin==3 %(r, 'tol', tol)
        tol = varargin{end}; 
    else %(r, sample_length)
        res = varargin{1}; 
    end       
end

poles = r.poles;  
dom = r.domain; 
a = dom(1); b = dom(2); 
%%
% the exponents are known in closed form via the poles of r.  
% We average the conjugate pairs of poles since conjugate symmetry is not 
% enforced by rfun:

%move poles to [0, 1) domain:
poles = (poles-a)/(b-a);
p = exp(-2i*pi*poles); 
p1 = p(abs(p)<1); 
p2 = p(abs(p)>1); 
p1 = sort(p1, 'ComparisonMethod', 'real');
p2 = 1./conj(p2); 
p2 = sort(p2, 'ComparisonMethod', 'real');

% average the poles: 
p = (p1 + p2)/2;

%%
% get a sample in Fourier space to solve for weights:
if isempty(res)
    res = r.res;
end
const = r.const; %normalize sample
scl = r.scl; 
r.const = 0; 
r.scl = 1; 
happy = false; 
while ~happy
    x = linspace(a, b, res+1); 
    x = x(1:end-1);
    y = sample2coeffs(feval(r,x), 'pos'); %does an fft
    %check for resolution: (check decay on coeffs):
    cutoff = standardChop(y, min(1e-10,tol));
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
        %s.nodes = log(p);  
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
        happy = norm(V*w-yy(xx+1)) < tol;
    else
        happy =true; 
    end
end

%%
% construct efun
 
s = efun(); 
s.exp = log(p);
s.weights = w; 
s.domain = r.domain; 
s.space = 'Fourier'; 
s.domain = r.domain; 
s.const = const; 
s.scl = scl;
res = res/2; 

s.res = ((res+mod(res,2))/2); 

end






