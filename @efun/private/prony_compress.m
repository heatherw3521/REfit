function [r, ss, w] = prony_compress(cp,m,tol)
%
% a version of Prony's method with rectangular Hankel matrices.

%%
chop_on = 1;

% length = 1 case: 
if ( length(cp)==1 || m==1 ) 
    %if the degree = 1, then the function should be the zero function:
    % the poles are arbitrary in this case. 
        w = 0; 
        r = log(.5);
        ss = [];
        return
end

if chop_on
    choptol = eps; 
    cutoff = standardChop(cp, choptol);
    n = max(cutoff, 2*m+1); %need at least 2m + 1 samples. 
    cp = cp(1:n); 
end

% build Hankel matrix  with m + 1 columns
n = length(cp); 
H = hankel(cp(1:n-m-1), cp(n-m-1:end)); 
[~, S, V]  = svd(H); 
ss = diag(S);
s = ss(ss > ss(1)*tol); 
M = min(m+1, length(s)+1);  
x = V(:,M);        
%%
% now find roots of Prony's polynomial: 
z = roots(flip(x(:).')); %solve via companion matrix
%get rid of roots at infinity and zero;
z = z(~isinf(z)); 
z = z(~(z==0));
%keep only roots that are in unit disk. 
z = z(abs(z) < 1); 

%% 
%roots of Prony polynomial are exponents for our sum. 
r = log(z); 
M = length(r);
%% find weights w/ LSQ sol to Vandermonde system

x = (0:n-1).';

V = repmat((z.'), n,1); 
V = V.^(repmat( x, 1,M));  
w = V\cp;
 
% shorten the sum if possible
if ~isempty(find(abs(w)<= tol, 1))
    idx = find(abs(w)> tol); 
    z = z(idx); 
    r = r(idx); 
    M = length(r); 
    %need to recompute the weights 
    V = repmat((z.'), n,1); 
    V = V.^(repmat(x, 1,M));  
    w = V\cp;
end

end


