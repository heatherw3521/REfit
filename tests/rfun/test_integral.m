function pass = test_integral()
%test rfun cumsum, integral, sum
tol = 1e-7; 
j = 1; 

f = @(x) exp(sin(2*pi*x))-1; 
x = linspace(0, 1, 2*500+2).'; x = x(1:end-1);
df = @(x) 2*pi.*cos(2*pi*x).*exp(sin(2*pi*x)); 

r = rfun(f(x), 'tol', 1e-8); 
dr = rfun(df(x), 'tol', 1e-8); 
r2 = cumsum(dr); 
pass(j) = max(abs(r2(x)-f(x))) < tol;
j = j+1; 

%test types:
r2 = cumsum(dr, 'rfun');
pass(j) = max(abs(r2(x)-f(x))) < tol;
j = j+1;

[r2, h] = cumsum(dr, 'rfun');
pass(j) = max(abs(r2(x)-f(x))) < tol;
j = j+1;
pass(j) = max(abs(h(x)-f(x))) < tol;
j = j+1;
pass(j) = isa(h, 'function_handle');
j = j+1; 

s2 = cumsum(dr, 'efun'); 
pass(j) = max(abs(s2(x, 'values')-f(x))) < tol;
j = j+1;

[s2, h]= cumsum(dr, 'efun');
pass(j) = isa(h, 'function_handle') && isa(s2, 'efun');
j = j+1;

% does function handle work when dr.const neq 0?
dr.const = 2; 
ff = @(x) f(x) + 2*x; 
h = cumsum(dr); 
pass(j) = max(abs(h(x)-ff(x))) < tol;
j = j+1;

% %test nonstandard domain: 
 f = @(x) exp(sin(pi*x))-1; %make it a mean zero function
 x = linspace(2, 4, 2*500+2).'; x = x(1:end-1);
 df = @(x) pi.*cos(pi*x).*exp(sin(pi*x)); 
 dr = rfun(df(x), 'tol', 1e-8, 'domain', [2,4]);
% 
 [r2, h] = cumsum(dr, 'rfun'); 
 pass(j) = max(abs(r2(x)-f(x))) < tol;
 j = j+1; 
 pass(j) = max(abs(f(x)-h(x))) < tol;
 j = j+1; 
% 
%check sum:
 pass(j) = abs(sum(dr))<1e-15; 
 j = j+1;
% 
 a= 2.5; b= 3.1; 
 I = sum(dr, [a, b]);
 pass(j) = abs(I -(f(b)-f(a)))<5*tol ; 
 j = j+1; 
% 
 pass(j) = abs(sum(dr, [2, 4]))<1e-15;
 j = j+1;
 pass(j) = abs(sum(dr, [2, 2]))==0;
% 
%%make sure integral works
 I = integral(dr, [a, b]);
 pass(j) = abs(I -(f(b)-f(a)))<5*tol ; 
 j = j+1; 
 
 dr.const = 2; %sum over total integral with non-zero const
 I = integral(dr); 
 pass(j) = abs(I -2*(diff(dr.domain)))<tol; 
 j = j+1;
% 
% %mean value zero function + const
 f = @(x) exp(sin(pi*x))-1+2*x;
 I = integral(dr, [a, b]); 
 pass(j) = abs(I -(f(b)-f(a)))<tol ;

end









