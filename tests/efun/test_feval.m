function pass = test_feval()
%test feval function
pass = []; 
% try on [0, 1] domain: 

fa = @(x) exp(sin(2*pi*(x+1)));
ss = efun(fa); 
tol = 1e-9; 
%vector of points in value space:
n = 1000; 
x = rand(n,1); 
j = 1;
pass(j) = max(abs(ss(x, 'values')-fa(x))) < tol;  % j = 1
j = j+1; 
pass(j) = max(abs(ss(x.', 'values')-fa(x))) < tol; %j = 2
j = j+1; 

% integers in frequency space: 
x = linspace(0,1, 1002); x = x(1:end-1).';
cf = sample2coeffs(fa(x)); 
l = length(cf); 
K = (-((l-1)/2):((l-1)/2)).';
pass(j) = max(abs(ss(K)-cf)) < tol;  %j = 3
j = j+1; 
pass(j) = max(abs(ss(K.')-cf)) < tol; %j = 4
j = j+1; 

% different domain (on [-1, 1]): 
fa2 = @(x) exp(sin(pi*(2*x-1))); 
x = linspace(-1, 1, 1000).';  
s = efun(fa2, 'domain', [-1, 1]); 
pass(j) = max(abs(s(x, 'values')-fa2(x))) < tol; %j = 5
j = j+1; 
pass(j) = max(abs(s(x.', 'values')-fa2(x))) < tol; %j = 6

end