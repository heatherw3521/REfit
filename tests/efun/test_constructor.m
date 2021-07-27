function pass = test_constructor()
%test basic functionality of constructor: 
j = 1; 
tol = 1e-8; 
%
%%test set 1: use spline function to test inputs/parsing: 
f1 = @(x) 3/4*abs(x).^3 - 3/2*(x).^2 + 1; 
f2 = @(x) 1/4*(2 - abs(x)).^3; 
f3 = @(x) f1(x).*(abs(x) <= 1) + f2(x).*(abs(x) > 1 & abs(x) <= 2); 
fa = @(x) f3(6*x-3); 
n = 500; xx = linspace(0, 1, 2*n+2); x = xx(1:end-1).'; 
%%
%input vals: j = 1
S = efun(fa(x)); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1; 
%input vals and locs: j = 2
S = efun(fa(x), x); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1; 
%input with bad sample set:  j= 3
t = rand(20,1); 
pass(j) = 0;
try
    S = efun(fa(t), t);
catch 
    pass(j) = 1;
end
j = j+1; 


%input a function handle (auto sample): j = 4
S = efun(fa); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;
%input a function handle with locs:  j = 5
S = efun(fa, x); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;
%input a function handle with sample j = 6: 
S = efun(fa, 2*n+1); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;
%input with sample # even: j = 7
S = efun(fa, 2*n+2); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;
%input with bad sample set: j = 8
pass(j) = 0;
try
    efun(fa, t);
catch 
    pass(j) = 1;
end
j = j+1; 

%input Fourier coeffs: j = 9
c = sample2coeffs(fa(x)); 
S = efun(c, 'coeffs'); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;

%pos only: j = 10
c = sample2coeffs(fa(x), 'pos'); 
xc = (0:(length(c)-1));
S = efun(c, xc, 'coeffs'); 
pass(j) = max(abs(S(x, 'values')-fa(x))) < tol;  
j = j+1;


%try with non-consecutive integers: j = 11
idx = sort(randperm(length(c), 301)).'; 
pass(j) = 0; 
try
    efun(c(idx), xc(idx), 'tol', 1e-8, 'coeffs'); 
catch
    pass(j) = 1; 
end
j = j + 1; 

% a few test with degree or tol specified: j = 12
m = length(S); 
pass(j) = 0;
try
    S = efun(fa(x), x, 'deg', m); 
catch 
    pass(j) = 1;
end
j = j+1; 

S = efun(fa,'tol', 1e-4); % j = 13
pass(j) = max(abs(S(x, 'values')-fa(x))) < 1e-3;  
j = j+1;

S = efun(fa(x),'tol', 1e-4); % j = 14
pass(j) = max(abs(S(x, 'values')-fa(x))) < 1e-3;  
j = j+1;

S = efun(c,(0:(length(c)-1)), 'tol', 1e-4, 'coeffs'); % j = 15
pass(j) = max(abs(S(x, 'values')-fa(x))) < 1e-3;  
j = j+1;

% try on a different domain:
fa = @(x) exp(sin(pi*(x))); 
dom = [-1, 1]; 
x = linspace(-1, 1, 1002).'; x = x(1:end-1); 
s = efun(fa, 'domain', dom); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < tol; 
j = j+1; 

s = efun(fa, 'domain', dom, 'tol', 1e-4); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < 1e-3; 
j = j+1; 

s = efun(fa(x), 'domain', dom); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < tol; 
j = j+1;

s = efun(fa(x), 'domain', dom, 'tol', 1e-4); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < 1e-3; 
j = j+1; 

s = efun(fa(x), x, 'domain', dom, 'tol', 1e-4); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < 1e-3; 
j = j+1;

s = efun(fa, 501, x, 'domain', dom, 'tol', 1e-4); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < 1e-3; 
j = j+1;

%test passing an rfun to constructor:
r = rfun(fa(x), x, 'tol', 1e-10, 'domain', dom); 
s = efun(r); 
pass(j) = max(abs(s(x, 'values')-fa(x))) < tol; 
j = j+1; 

%test passing an efun for compression: 
s2 = efun(s,'tol', 1e-4); 
pass(j) = max(abs(s2(x, 'values')-fa(x))) < 1e-4; 
j = j+1; 
pass(j) = length(s2) <= length(s); 

end







