function s = ft(r, varargin)
% Compute the Fourier transform of an rfun. 
% This returns an efun object (sum of complex exponentials)
% representing the Fourier transform of the trigonometric rational r.
%
%%
% See also: efun/efun, efun/ift
%
% For more information, see Wilber, H., Damle, A., and Townsend, A. (2021). 



%%
%
if isempty(r)
    s = [];
    return
end

if all(r.vals==0) %deal with zero function and constants
        s = efun([0;0]);
        s.domain = r.domain;
        s.const = r.const;
        s.scl = r.scl;
        return
end
    
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
    res = max(r.res+(1-mod(r.res, 2)), 251);
end
happy = false; 
while ~happy
    x = linspace(a, b, res+1); 
    x = x(1:end-1);
    y = sample2coeffs(feval(r,x), 'pos'); %does an fft
    %check for resolution: (check decay on coeffs):
    cutoff = standardChop(y, min(1e-10,tol));
    happy = cutoff < length(y); 
    if res > 20000 % max sample
        break
    end
    res = 2*res; 
end
%Ns = length(y)-1;
%yy = y(1:max(cutoff, round(Ns/2))); %sometimes cutoff can be too severe, this is a bandaid for now
yy = y(1:cutoff); 
N = length(yy);  
const = yy(1); 
yy(1) = 0; 
scl = max(abs(yy)); 
yy = yy/scl; 
%%
% use a subset of the samples to solve for weights:
mm = length(poles); 
happy = false; 
while ~happy 
    mm = min(2*mm, N); 
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
    sp = min(N, 20); 
    xx = [(0:min(50, mm-1)).'; randi([0, N-1], sp,1)]; 
    nn = length(xx); 
    V = repmat((p.'), nn,1); 
    V = V.^(repmat(xx, 1,length(p)));  
    happy = norm(V*w-yy(xx+1)) < tol;
    if ~happy && mm == N
        %if ft fails give up and call efun on the sample:
        yy = scl*yy;
        yy(1) = const; 
        yy = [flip(conj(yy(2:end)));yy]; 
        s = efun(yy, -(N-1):(N-1), 'tol', tol, 'domain', dom, 'coeffs');
        return 
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
s.tol = r.tol;
 

s.res = ((N+mod(N,2))/2); 

end






