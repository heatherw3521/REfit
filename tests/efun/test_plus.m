function pass = test_plus()
%test plus function
pass = []; 
[s, fa] = REFIT.gallerysum('chirp');
x = linspace(0, 1, 2*1000 + 2); 
x = x(1:end-1).'; 
tol = 1e-8; 
j = 1;
% plus self
s1 = 2*s; 
s2 = s+s; 
pass(j) = max(abs(s1(x, 'values')-2*fa(x))) < 1.5*tol;  % j = 1
j = j+1; 
pass(j) = max(abs(s2(x, 'values')-s1(x, 'values'))) < 1e-15;  % j = 2
j = j+1; 

% plus constant: 
s1 = s+10; 
s2 = efun(fa(x)+10, 'tol', 1e-9); 
pass(j) = max(abs(s2(x, 'values')-s1(x, 'values'))) < tol;  
j = j+1; 

%plus two efuns: 
[g, fb] = REFIT.gallerysum('wild');
%h = efun(fb(x) + fa(x), 'tol', 1e-10);
h = @(x) fb(x) + fa(x); 
h2 = s + g; 
pass(j) = max(abs(h2(x, 'values')-h(x))) < 3*tol;  

%plus an rfun and efun:



end