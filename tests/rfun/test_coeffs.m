function pass = test_coeffs()
% test the coeffs function
tol = 1e-9; 
j = 1; 
pass(j) = isempty(coeffs(rfun())); 
j = j + 1; 

[r, fa] = gallery_rfun('wild'); 
p = r.res; 
x = linspace(0, 1, p+1); x = x(1:end-1)'; 
cf = sample2coeffs(fa(x), 'pos'); 
cr = coeffs(r); 
n = length(cr); 
m = (n-1)/2+1;
cf = cf(1:m); 
cf = [flip(conj(cf(2:end))); cf];

pass(j) = max(abs(cr-cf)) < 1e-6;
j = j + 1; 

modes = (-(m-1):(m-1)).'; 
pass(j) = max(abs(coeffs(r, modes)-cf)) < 1e-6;

end

