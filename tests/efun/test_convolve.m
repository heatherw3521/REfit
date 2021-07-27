function pass = test_convolve()
%test convolve
tol = 1e-8; 
j = 1; 

[s1, f1] = gallery_efun('spline'); 
%convolve with self as a simple test: 

tru = @(x) f1(x).*f1(x); 
x = linspace(0, 1, 500).'; x = x(1:end-1); 
s3 = convolve(s1, s1); 
pass(j) = max(abs(s3(x, 'values') - tru(x)))<tol; 
j = j+1; 

%try an rfun as 2nd input:
r1 = ift(s1);
s3 = convolve(s1, r1);
pass(j) = max(abs(s3(x, 'values') - tru(x)))<tol; 
j = j+1; 

 
%try function handle as 2nd input:
f = @(n) s1(n); 
s3 = convolve(s1, f);
pass(j) = max(abs(s3(x, 'values') - tru(x)))<tol; 
j = j+1;

%try on conflicting domains: 
s2 = s1; 
s2.domain = [3, 4]; 
pass(j) = 0;
try
    s3 = convolve(s1, s2);
catch 
    pass(j) = 1;
end

% test a convolution that requires more poles
[r2, f2] = gallery_rfun('wild'); 
s2 = ft(r2);
const = s2.const;
s2.const = 0;
r2.const = 0; 
f2 = @(x) f2(x) - const; 

tru = @(x) f1(x).*f2(x);  
s3 = convolve(s1, s2); 
pass(j) = max(abs(s3(x, 'values') - tru(x)))<tol; 

end



