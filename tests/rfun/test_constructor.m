function pass = test_constructor()
% test the rfun constructor: 
warning off
f1 = @(x) abs(sin(2*pi*x)); 
f2 = @(x) abs(sin(2*pi*x)).^4; 
f3 = @(x) abs(sin(pi*x)); 
f4 = @(x) abs(sin(pi*x)).^4; 
tol = 1e-9;

%NOTE: f2/f4 are for testing the auto-sample feature.
% For stronger singularities, one needs to sample on a grid
% that clusters around the singularities. The autosampling is 
% equally spaced and doesn't have good convergence properties, so 
% it isn't an especially useful feature right now. Future work could 
% include identifying singularity locations as AAA is executed 
% and adaptively refining the grid appropriately. 

%try on standard domain:
pass1 = test_it([0, 1],f1, f2, tol);

%try on domain [2, 4]:
pass2 = test_it([2, 4], f3, f4, tol); 

pass = [pass1, pass2]; 
j = length(pass) + 1; 

%check some basic things about r properties: 
n = 2000; 
tol = 1e-8; 
x = linspace(3, 5, 2*n+2).'; x = x(1:end-1); 
r = rfun(f3(x),x, 'dom', [3, 5], 'tol', tol); 
pol = r.poles; 
pass(j) = (all(real(pol) >=3) && all(real(pol) < 5)); %pass(33)
j = j+1; 
pass(j) = (length(pol)==length(r)); %34
j = j+1;
pass(j) = (length(r.nodes)==length(r)); %35
j = j+1; 
pass(j) = (length(r.weights) == length(r)); %36

warning on
end

%%%%%%%%%
function pass = test_it(dom, f1, f2, tol)

a = dom(1); 
b = dom(2); 
if all(dom == [0, 1])
    dom = []; %test the default
end

n = 2000; x = linspace(a, b, 2*n+2); x = x(1:end-1).';
j = 1;

% test vector:
if isempty(dom)
    r = rfun(f1(x), x);
else
    r = rfun(f1(x), x, 'dom', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %pass(1)
j = j+1; 

if isempty(dom)
    r = rfun(f1(x), x, 'tol', tol); 
else 
    r = rfun(f1(x), x, 'tol', tol, 'domain', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %2
j = j+1; 

deg = length(r); 
if isempty(dom)
    r = rfun(f1(x), x, 'deg',deg); 
else
    r = rfun(f1(x), x, 'deg',deg, 'dom', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %3
j = j+1; 

if isempty(dom)
    r = rfun(f1(x), x, 'mmax',deg); 
else
    r = rfun(f1(x), x, 'mmax',deg, 'dom', dom);
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %4
j = j+1;

if isempty(dom)
    r = rfun(f1(x), x, 'tol', 1e-4, 'mmax',deg); 
else
    r = rfun(f1(x), x, 'dom', dom, 'tol', 1e-4, 'mmax',deg); 
end
pass(j) = max(abs(r(x)-f1(x))) < 1e-3; %5
j = j+1;

if isempty(dom)
    r = rfun(f1(x), x, 'tol', 1e-4, 'mmax',deg, 'cleanup_off'); 
else
    r = rfun(f1(x), x, 'tol', 1e-4, 'dom', dom, 'mmax',deg, 'cleanup_off');
end
pass(j) = max(abs(r(x)-f1(x))) < 1e-3; %6
j = j+1;

% test vector only: 
if isempty(dom)
    r = rfun(f1(x), 'tol', tol); 
else
    r = rfun(f1(x), 'domain', dom, 'tol', tol);
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %7
j = j+1; 

if isempty(dom)
    r = rfun(f1(x), 'deg', deg); 
else
    r = rfun(f1(x), 'deg', deg, 'dom', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %8
j = j+1; 

if isempty(dom)
    r = rfun(f1(x), 'deg', deg, 'tol', 1e-4); 
else
    r = rfun(f1(x), 'deg', deg, 'dom', dom, 'tol', 1e-4);
end
pass(j) = max(abs(r(x)-f1(x))) < 1e-3; %9
j = j+1; 

if isempty(dom)
    r = rfun(f1(x), 'cleanup_off'); 
else
    r = rfun(f1(x), 'cleanup_off', 'domain', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < 1e-3; %10
j = j+1; 


%test function handle: 
if isempty(dom)
    r = rfun(f1, x); 
else 
    r = rfun(f1, x, 'dom', dom); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %11
j = j+1; 

if isempty(dom)
    r = rfun(f1, x, 'tol', tol); 
else
    r = rfun(f1, x, 'domain', dom, 'tol', tol); 
end
pass(j) = max(abs(r(x)-f1(x))) < tol; %12
j = j+1; 


if isempty(dom)
    r = rfun(f1, 2000); %input sample length
    [~, s] = sample(r); %get the sample back
else
    r = rfun(f1, 2000, 'dom', dom);
    [~, s] = sample(r);
end
pass(j) = max(abs(r(s)-f1(s))) < 2*tol; %13
j = j+1;

if isempty(dom)
    r = rfun(f1, 2000, 'mmax', deg); 
    [~, s] = sample(r); %get the sample back
else
    r = rfun(f1, 2000, 'mmax', deg, 'domain', dom); 
    [~, s] = sample(r); %get the sample back
end
pass(j) = max(abs(r(s)-f1(s))) < tol; %14
j = j+1;

%test auto-sample: 
if isempty(dom)
    r = rfun(f2); 
else
    r = rfun(f2, 'dom', dom); 
end
pass(j) = max(abs(r(x)-f2(x))) < 1e2*tol; %15
j = j+1; 

%try a different tolerance
if isempty(dom)
    r = rfun(f2, 'deg', deg, 'tol', 1e-4); 
else
    r = rfun(f2, 'domain', dom, 'deg', deg, 'tol', 1e-4); 
end
pass(j) = max(abs(r(x)-f2(x))) < 1e-3; %16

end

