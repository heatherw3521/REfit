function pass = test_poles()
%test basic functionality  of get.poles

[f, s] = gallery_efun('spline'); 
x = linspace(0, 1, 3002).'; x = x(1:end-1); 
tol = 1e-8; 
j = 1; 
%1: basic construction: 
r = rfun(f(x), x, 'tol', tol); 
p = r.poles; 
pass(j) = (length(p)==length(r)); 
j = j+1; 

pol = poles(r); 
pass(j) = all((pol-p)==0); 
j = j+1;

%poles should be off real line: 
pass(j) = all(abs(imag(pol)) > 1e-10); 
j = j + 1; 
%poles should be in [0, 1]:
rp = real(r.poles); 
pass(j) = all( rp >= 0 & rp <= 1); 
j = j+1;

pol = poles(r, 'zt'); 
pass(j) = all(abs(pol)<1); 
j = j+1; 
pass(j) = (length(pol) == length(r)/2); 
j = j+1; 
%try alternative domain 
r.domain = [4, 5]; 
rp = real(r.poles); 
pass(j) = all( rp >= 4 & rp <= 5); 
j = j+1;

% now test when constructed via inverse Fourier transform: 
r2 = ift(s); 
pass(j) = (length(r2.poles)==2*length(s));
j = j+1; 
pol = poles(r2, 'zt'); 
pass(j) = (length(pol) == length(s)); 
j = j+1; 
pass(j) = all(abs(pol)< 1); 
j = j+1; 
pass(j) = all(abs(imag(r2.poles)) > 1e-10);
rp = real(r2.poles); 
pass(j) = all(rp >= 0 & rp <= 1); 
j = j+1;
%test different domain for s: 
s.domain = [4, 5]; 
r2= ift(s); 
pass(j) = (length(r2.poles)==2*length(s));
j = j+1; 
rp = real(r2.poles); 
pass(j) = all(rp >= 4 & rp < 5); 

end




