function pass = test_diff()
% test efun diff function: basic test of functionality
tol = 1e-9; 
j = 1; 

f = @(x) exp(sin(2*pi*x))-1.266065877752008; %make it a mean zero function
x = linspace(0, 1, 2*500+2).'; x = x(1:end-1);
df = @(x) 2*pi.*cos(2*pi*x).*exp(sin(2*pi*x)); 
df2 = @(x) -4*pi^2.*sin(2*pi*x).*exp(sin(2*pi*x)) + ...
    (2*pi.*cos(2*pi*x)).^2.*exp(sin(2*pi*x));

s = efun(f(x), 'tol', 1e-11); 
s2 = efun(df(x), 'tol', 1e-11); 
s3 = efun(df2(x), 'tol', 1e-11); 
ds = diff(s); 
idx = (0:20).'; 
pass(j) = max(abs(s2(idx)-ds(idx))) < tol; 
j = j+1; 

ds = diff(s, 'efun'); 
pass(j) = max(abs(df(x)-ds(x, 'values'))) < tol;
j = j+1;

%rfun + handle
[rsd, hh] = diff(s, 'rfun'); 
pass(j) = max(abs(rsd(x)-df(x))) < 5e1*tol; %lose some accuracy here.
j = j+1;
pass(j) = max(abs(ds(idx)-hh(idx))) < tol; %lose some accuracy here.
j = j+1;

dss = diff(s, 2); 
pass(j) = max(abs(s3(idx)-dss(idx))) < tol; 
j = j+1;

%efun + handle
[dss, dsh] = diff(s, 'efun', 2); 
pass(j) = max(abs(df2(x)-dss(x, 'values'))) < tol;
j = j+1;
pass(j) = max(abs(s3(idx)-dsh(idx))) < tol; 
j = j+1; 

dss = diff(s, 2, 'rfun'); 
pass(j) = max(abs(df2(x)-dss(x))) < tol;

%try on a nonstandard domain: 
f = @(x) exp(sin(pi*x))-1.266065877752008; %make it a mean zero function
x = linspace(2, 4, 2*500+2).'; x = x(1:end-1);
df = @(x) pi.*cos(pi*x).*exp(sin(pi*x)); 
s =  efun(f(x), 'tol', 1e-11, 'domain', [2,4]);
s2 = efun(df(x), 'tol', 1e-11, 'domain', [2,4]); 

[ds, h] = diff(s, 'efun'); 
pass(j) = max(abs(s2(idx)-h(idx))) < tol; 
j = j+1;
pass(j) = max(abs(df(x)-ds(x, 'values'))) < tol;
j = j+1;

%%
% try an example with slower F.C. decay and a looser tol.
% Here the F.C. of the derivative decay much slower than those of f. 
% this dominates the error. 
tol = 1e-3; 
x = linspace(0, 1, 2*3000+2).'; x = x(1:end-1); 
f = @(x) abs(sin(2*pi*x)).^3-0.424342469618778; 
df = @(x)  3*pi*sin(4*pi*x).*abs(sin(2*pi*x)); 
s = efun(f(x), 'tol', 1e-5); 
s2 = efun(df(x), 'tol', 1e-5); 
[ds, h] = diff(s, 'rfun');
idx = (0:100).'; 
pass(j) = max(abs(s2(idx)-h(idx))) < tol; 
j = j+1; 
pass(j) = max(abs(df(x)-ds(x, 'values'))) < 5*tol; 

end








