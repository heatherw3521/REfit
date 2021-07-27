function pass = test_convolve()
%test convolve
tol = 1e-8; 
j = 1; 

[s1, f1] = gallery_efun('spline'); 
%convolve with self as a simple test: 
r = ift(s1); 
tru = s1.*s1; 
x = linspace(0, 1, 500).'; x = x(1:end-1); 
r2 = convolve(r, r); 
pass(j) = max(abs(r2(x) - tru(x, 'values')))<tol; 
j = j+1; 

%try an efun as 2nd input:
r2 = convolve(r, s1);
pass(j) = max(abs(r2(x) - tru(x, 'values')))<tol; 
j = j+1; 


%try function handle as 2nd input:
r2 = convolve(r, f1);
pass(j) = max(abs(r2(x) - tru(x, 'values')))<tol; 
j = j+1;

%try on conflicting domains: 
rr = r; 
rr.domain = [3, 4]; 
pass(j) = 0;
try
    r2 = convolve(r, rr);
catch 
    pass(j) = 1;
end

% test a convolution that requires more poles
[r2, f2] = gallery_rfun('wild'); 
const = r2.const;
r2.const = 0;
f2 = @(x) f2(x) - const; 
s2 = ft(r2); 

tru = s2.*s1;   
r3 = convolve(r, r2); 
pass(j) = max(abs(r3(x) - tru(x, 'values')))<tol; 

end



