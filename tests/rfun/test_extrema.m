function pass = test_extrema()
% test basic functionality of min, max. 

tol = 1e-4; 
j = 1; 

f = @(x) abs(sin(2*pi*x)).^3; 
x = linspace(0, 1, 8000).'; x = x(1:end-1); 
r = rfun(f(x), x, 'tol', 1e-5); 

tmx = [1/4,3/4].'; 
tmn = [1/2,1].'; 

[mn, mx] = extrema(r); 
pass(j) = max(abs(sort(mx)-tmx)) < tol; 
j = j+1; 
pass(j) = max(abs(sort(mn)-tmn))<2*tol;
j = j+1; 

%try on different domain: 
r.domain = [2, 3]; 
tmx = tmx+2; 
tmn = tmn+2; 
[mn, mx] = extrema(r); 
pass(j) = max(abs(sort(mx)-tmx)) < tol; 
j = j+1; 
pass(j) = max(abs(sort(mn)-tmn))<2*tol;
j = j+1; 

[s, fa] = gallery_efun('spline'); 
s2 = s.*(-s); %convolve the rationals assoc. with s, -s.
r = ift(s2); 
mx = max(r, 'abs'); 
pass(j) = abs(mx-.5)< tol; 
mn = min(-r, 'abs'); 
pass(j) = abs(mn-.5)< tol;


end


