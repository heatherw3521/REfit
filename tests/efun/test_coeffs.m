function pass = test_coeffs()
% test the getcoeffs function

tol = 1e-9; 
j = 1; 
pass(j) = isempty(coeffs(efun())); 
j = j + 1; 

[s, fa] = REFIT.gallerysum('wild'); 
p = s.res; 
x = linspace(0, 1, 2*p + 2); x = x(1:end-1)'; 
cf = sample2coeffs(fa(x)); 

pass(j) = max(abs(coeffs(s)-cf)) < tol;
j = j + 1; 

modes = (-p:p).'; 
pass(j) = max(abs(coeffs(s, modes)-cf)) < tol;

end

