function pass = test_integral()
%test efun cumsum, integral, sum
tol = 1e-9; 
j = 1; 

f = @(x) exp(sin(2*pi*x))-1.266065877752008; %make it a mean zero function
x = linspace(0, 1, 2*500+2).'; x = x(1:end-1);
df = @(x) 2*pi.*cos(2*pi*x).*exp(sin(2*pi*x)); 
df2 = @(x) -4*pi^2.*sin(2*pi*x).*exp(sin(2*pi*x)) + ...
    (2*pi.*cos(2*pi*x)).^2.*exp(sin(2*pi*x));

s = efun(f(x), 'tol', 1e-8); 
ds = efun(df(x), 'tol', 1e-8); 
s2 = cumsum(ds); 
idx = (0:100).'; 
pass(j) = max(abs(s(idx)-s2(idx))) < tol;
j = j+1; 

%test outputs: 
[s2, h] = cumsum(ds, 'efun'); 
pass(j) = max(abs(s(idx)-h(idx))) < tol;
j = j+1; 
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < tol;
j = j+1; 
s2 = cumsum(ds, 'efun'); 
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < tol;
j = j+1; 
r2 = cumsum(ds, 'rfun'); 
pass(j) = max(abs(s(x, 'values')-r2(x))) < tol;

%test nonstandard domain: 
f = @(x) exp(sin(pi*x))-1.266065877752008; %make it a mean zero function
x = linspace(2, 4, 2*500+2).'; x = x(1:end-1);
df = @(x) pi.*cos(pi*x).*exp(sin(pi*x)); 
s =  efun(f(x), 'tol', 1e-8, 'domain', [2,4]);
ds = efun(df(x), 'tol', 1e-8, 'domain', [2,4]);

[s2, h] = cumsum(ds, 'efun'); 
pass(j) = max(abs(s(idx)-h(idx))) < tol;
j = j+1; 
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < tol;
j = j+1; 

%check sum:
pass(j) = abs(sum(ds))<1e-15; 
j = j+1;

a= 2.5; b= 3.1; 
I = sum(ds, [a, b]);
pass(j) = abs(I -(s(b,'values')-s(a, 'values')))<5*tol ; 
j = j+1; 

pass(j) = abs(sum(ds, [2, 4]))<1e-15;
j = j+1;
pass(j) = abs(sum(ds, [2, 2]))==0;

%make sure integral works
I = integral(ds, [a, b]);
pass(j) = abs(I -(s(b,'values')-s(a, 'values')))<5*tol ; 
j = j+1; 

ds.const = 1; %sum over total integral with non-zero const
I = integral(ds); 
pass(j) = abs(I -2)<tol; 
j = j+1;

%mean value zero function + const
f = @(x) exp(sin(pi*x))-1.266065877752008+x;
I = integral(ds, [a, b]); 
pass(j) = abs(I -(f(b)-f(a)))<tol ;

end









