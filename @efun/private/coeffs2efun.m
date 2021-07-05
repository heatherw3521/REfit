
function [w,r, res, ss]= coeffs2efun(cp, x, chop_on, tol)
%%
%  cp is a vector of Fourier coefficients, x are the modes (modes are geq 0). 
%  Returns mx1 vectors of weights (w) and nodes (r), 
%  so that cp_(x(k)) \approx \sum_{j = 1}^{m}w_j exp(r(x(k))) 
%
%  coeffs2exp(f, tol, 'rand') computes w and r using a randomized SVD,
%  which is faster. 
%
%  If deg is populated, then the exponential sum is built so that m \leq
%  deg. 

%%
cp = cp(:); 
n = length(cp); 

%%
% it can hurt the accuracy of the method if we 
% are trying to fit a lot of coeff values that are 
% small enough to be considered noise. 
%
% we borrow chebfun's 'chop' function to truncate
% the coeffs: 

if chop_on
choptol = eps; 
cutoff = standardChop(cp, choptol);
 
%make the cutoff have the right parity:
N = cutoff + mod(cutoff+1, 2);
N = N + 2*mod((N-1)/2, 2); 
if N > n %pad with zeros
    cp = [cp; zeros(N-n,1)]; 
    x = [x; ((x(end)+1):x(end)+N-n).'];
else
    cp = cp(1:N); 
    x = x(1:N); 
end
%%
% find roots of Prony polynomial:
if ( length(cp)==1 )
    %if the degree = 1, then the function should be the zero function:
    % the poles are arbitrary in this case. 
        w = 0; 
        r = log(.5);
        res = 1;
        ss = [];
        return
end
[z,ss] = coneigen(cp,tol); 
end

%roots of Prony polynomial are exponents for our sum. 
r = log(z); 
M = length(r);
%% find weights w/ LSQ sol to Vandermonde system
if isempty(x)
    x = (0:N-1).';
end
V = repmat((z.'), N,1); 
V = V.^(repmat( x, 1,M));  
w = V\cp;
 
% shorten the sum if possible
if ~isempty(find(abs(w)<= tol, 1))
    idx = find(abs(w)> tol); 
    z = z(idx); 
    r = r(idx); 
    M = length(r); 
    %need to recompute the weights 
    V = repmat((z.'), N,1); 
    V = V.^(repmat(x, 1,M));  
    w = V\cp;
end
res = N-1; %return max mode. 
end 










    















