function [z, ss] = coneigen(cp,tol, svd_type)
%
% finds z, the roots of the Prony polynomial associated
% with the vector cp using a variation of Prony's method. 
%
% see Beylkin,G.,  Monzon, L. (2009)
% 
if isempty(svd_type)
    svd_type = 'randomized'; 
end

%%
% build a square Hankel matrix from these coeffs 
n = length(cp); %needs to be odd. 
m = (n+1)/2; %size of Hankel matrix (square)

 
switch svd_type
    case 'standard'
    H = hankel(cp(1:m), cp(m:end)); 
    [~, S, V]  = svd(H); 
    ss = diag(S);
    s = ss(ss > ss(1)*tol); 
    M = length(s);  
    if M == length(ss) 
        fprintf('No Hankel sing value is small enough, acc will be %e\n',s(end))
    else
        M = M+1; %pick vector in the approximate null space. 
    end
    x = V(:,M); 
        
    case 'randomized'
        cl = cp(1:m); 
        rw = cp(m:end); 
        p=10; %oversampling
        l = min(10, m-1-p);  %guess the rank 
        [x,ss] = coneigen_rand(cl,rw, l, tol);   
end

%%
% now find roots of Prony's polynomial: 
z = roots(flip(x(:).')); %solve via companion matrix
%get rid of roots at infinity and zero;
z = z(~isinf(z)); 
z = z(~(z==0));
%keep only roots that are in unit disk. 
z = z(abs(z) < 1); 
end

