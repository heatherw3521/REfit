function pass = test_roots()
j = 1; 
tol = 1e-8; 
f = @(x) abs(sin(2*pi*x)).^3-.5; 
zer = [asin((.5)^(1/3))/2/pi, asin(-(.5)^(1/3))/2/pi] ; 
zer = mod([zer,  zer + 1/2], 1).'; 
zer = sort(zer); 
x = linspace(0, 1, 8000).'; x = x(1:end-1); 

r = rfun(f(x), x, 'tol', 1e-8); 

rt = roots(r); 

pass(j) = abs(sort(rt)-sort(zer))< tol; 
j = j+1;

%check for all roots: 
rt = roots(r, 'all');
pass(j) = (length(rt) == length(r)); 
%now check with zero constant: should be type n-1, n:
r.const = 0; 
rt = roots(r, 'all'); 
pass(j) = (length(rt) == length(r)-2); 
j = j+1; 

% try different domain: 
f = @(x) abs(sin(pi*x)).^3-.5; 
zer = [asin((.5)^(1/3))/pi, asin(-(.5)^(1/3))/pi] ; 
zer = mod([zer,  zer + 1], 2).'; 
zer = zer + 2; 
zer = sort(zer); 
x = linspace(2, 4, 8000).'; x = x(1:end-1); 
r = rfun(f(x), x, 'tol', 1e-8, 'domain', [2, 4]); 
rt = roots(r); 

pass(j) = abs(sort(rt)-sort(zer))< tol; 
j = j+1;
%check for all roots: 
rt = roots(r, 'all');
pass(j) = (length(rt) == length(r)); 
%now check with zero constant: should be type n-1, n:
r.const = 0; 
rt = roots(r, 'all'); 
pass(j) = (length(rt) == length(r)-2); 


end