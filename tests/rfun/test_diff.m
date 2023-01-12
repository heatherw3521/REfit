function pass = test_diff()
% test rfun diff function: basic test of functionality
%
tol = 1e-8; 
j = 1; 

f = @(x) exp(sin(2*pi*x)); 
x = linspace(0, 1, 2*500+2).'; x = x(1:end-1);
df = @(x) 2*pi.*cos(2*pi*x).*exp(sin(2*pi*x)); 
df2 = @(x) -4*pi^2.*sin(2*pi*x).*exp(sin(2*pi*x)) + ...
    (2*pi.*cos(2*pi*x)).^2.*exp(sin(2*pi*x));
%%
r = rfun(f(x),x, 'tol', 1e-9); 
dr = diff_bary(r, 2); 
pass(j) = max(abs(dr(x)-df(x))) < 2*tol; 
j = j+1; 

dr = diff(r, 'efun'); 
pass(j) = max(abs(dr(x, 'values')-df(x))) < 5e1*tol;
j = j+1;

%rfun + handle
[dr, hh] = diff(r, 'rfun'); 
pass(j) = max(abs(dr(x)-df(x))) < 5e1*tol; %lose some accuracy here.
j = j+1;
pass(j) = max(abs(dr(x)-hh(x))) < 2*tol; %lose some accuracy here.
j = j+1;
pass(j) = isa(hh, 'function_handle');
j = j+1; 

%try higher derivative
dr2 = diff(r, 2); 
pass(j) = max(abs(dr2(x)-df2(x))) < 5e2*tol; 
j = j+1;

%efun + handle
[dr2, hh2] = diff(r, 'efun', 2); 
pass(j) = max(abs(df2(x)-dr2(x, 'values'))) < 2e2*tol;
j = j+1;
pass(j) = max(abs(hh2(x)-df2(x))) < 2e2*tol; 
j = j+1; 
pass(j) = isa(hh2, 'function_handle');
j = j+1; 

dr2 = diff(r, 2, 'rfun'); 
pass(j) = max(abs(df2(x)-dr2(x))) < 2e2*tol;
j = j+1;

%try constant function: 
r = rfun([0; 0], [1/4; 3/4]); 
r.const = 2; 
dr =  diff(r); 
pass(j) = all(dr(x)==0); 
j = j+1;
[dr, h] = diff(r, 'efun'); 
pass(j) = isa(h, 'function_handle'); 
j = j+1; 
pass(j) = all(dr(0:10)==0) && isa(dr, 'efun'); 
j = j+1; 
dr = diff(r, 3, 'rfun'); 
pass(j) = all(dr(x)==0)&& isa(dr, 'rfun'); 


%try on a nonstandard domain: 
f = @(x) exp(sin(pi*x)) ;
x = linspace(2, 4, 2*500+2).'; x = x(1:end-1);
df = @(x) pi.*cos(pi*x).*exp(sin(pi*x)); 
r =  rfun(f(x), 'tol', 1e-9, 'domain', [2,4]);

[dr, h] = diff(r, 'rfun'); 
pass(j) = max(abs(dr(x)-df(x))) < 5e1*tol; 
j = j+1;
pass(j) = max(abs(df(x)-h(x))) < 5*tol;
j = j+1;

%%
% try an example with slower F.C. decay and a looser tol.
% Once the closed-form formula is available, this should perform better.
tol = 1e-3; 
x = linspace(0, 1, 2*3000+2).'; x = x(1:end-1); 
f = @(x) abs(sin(2*pi*x)).^3-0.424342469618778; 
df = @(x)  3*pi*sin(4*pi*x).*abs(sin(2*pi*x)); 
r = rfun(f(x), 'tol', 1e-5);  
[dr, h] = diff(r, 'rfun'); 
pass(j) = max(abs(dr(x)-df(x))) < 5*tol; 
j = j+1; 
pass(j) = max(abs(df(x)-h(x))) < 5*tol; 

end








