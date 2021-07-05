
function [w,r, ss]= vals2efun(f,x, varargin)
%%
%  vals2efun(f, x,tol), where f(x) is a vector of values at locations x, 
%  returns m by 1 vectors of weights (w) and nodes (r), 
%  so that the values are well approximated by the sum f(v) \approx \sum_{j = 1}^w_j exp(r(v)) 
%

%% 
% determine tolerance
tol = 1e-8; 
[~, tolid] = find(strcmp(varargin,'tol'));
if ~isempty(tolid)
    tol = varargin{tolid+1}; 
end
 
f = f(:); 
n = length(f); 
x = x(:); 

%%
deg = []; 
[~, degid] = find(strcmp(varargin,'degree'));
if ~isempty(degid)
    deg = varargin{degid+1}; 
end

%if length is even, throw error: 
if mod(n, 2)==0
    error('esum:vals2exp:vector of values must be odd in length.')
end

if ~isempty(deg)
    [z,ss] = coneigen_APM_rect(f,deg);
elseif any(strcmp(varargin, 'rand'))
    l = min(10, max( floor((length(cp) - 3)/2)-10,1)); %we guess an initial rank/order of the exp sum.  
    z = coneigen_rand(f, l, tol); 
else 
    [z,ss] = coneigen_APM(f, tol); 
end

%roots of Prony polynomial are exponents for our sum. 
r = (n-1)*log(z); %the n-1 is for scaling the interval. won't work if points aren't eq. sp.
M = length(r); 
%r = log(z); 
%M = length(r); 
%% find weights w/ LSQ sol to Vandermonde system

V = repmat((z.'), n,1); 
V = V.^(repmat( x*(n-1), 1,M));  
%V = V.^(repmat(x, 1, M)); 
w = V\f;


% shorten the sum if possible
 if ~isempty(find(abs(w)<= tol, 1)) && isempty(degid)
     idx = find(abs(w)> tol); 
     z = z(idx); 
     r = r(idx); 
     M = length(r); 
     %need to recompute the weights 
     V = repmat((z.'), n,1); 
     V = V.^(repmat( x*(n-1), 1,M));  
     w = V\f;
 end 

%check for degree: %if needed, we add in zero weights: 
if ~isempty(degid) && ~(length(w) == deg)
    K = deg -length(w); 
    w = [w; zeros(K,1)];
    r = [r; -100*ones(K,1)]; 
end   

end 










    















