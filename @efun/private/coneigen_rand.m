function [x,s] = coneigen_rand(cl,rw, l, tol)
%
% Uses a randomized svd method (as in Halko, Martinsson, Tropp) 
% to find the svd(con-eigen pairs) of the Hankel matrix H 
% associated with f and then find roots of prony polynomial. 
%
% takes advantage of fast Hankel matvec. 

n1 = length(rw); 
%n2 = length(cl); 
p = 10; % oversampling parameter

while l+p < n1
    [s,v] = rand_svd_hankel(cl,rw,n1,l+p);
    % TO DO: fast update algorithm as l increases.
    % right now we recompute everything. 
    s = diag(s);
    s2 = abs(s)/max(abs(s));
    idx = find(s2<tol,1);
    if isempty(idx)
        l = l+10;
    else
        break;
    end
end
if isempty(idx)
	fprintf('No coneigen value is small enough, acc will be %e\n',s2(end))
    %x = v(:,end)+conj(u(:,end));
    x = v(:, idx);
else
    %x = v(:,idx)+conj(u(:,idx));
    x = v(:, idx); 
end 
end

function [S,V] = rand_svd_hankel(cl,rw,N,l)
% USV approx = hankel(cl, rw). 
% (Here, U = Q*W, where W are the right sing. vecs of B)

%appx the range of H
Omega = randn(N,l);
Y = hankelmult(rw, cl, Omega); 
[Q, ~] = qr(Y,0);

%%
% appx svd. 
B = (hankelmult(conj(rw), conj(cl), Q))';
[~, S, V] = svd(B,'econ');
%U = Q*U; %note: this is not always the Takagi factorization. 
end

function y=hankelmult(a,b,x)
%a = last row
%b = first col
a = a(:); 
b = b(:);    
[~, m]= size(x);
x = flip(x); 
n=length(b);
c = [a; 0; (b(1:end-1))];
F=fft(c);
t=ifft(F.*fft([x; zeros(n,m)]));
y = t(1:n,:);  
end
