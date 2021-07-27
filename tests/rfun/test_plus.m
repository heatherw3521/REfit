function pass = test_plus()
%test plus function 
[r, fa] = gallery_rfun('wells');
x = linspace(0, 1, 2*1000 + 2); 
x = x(1:end-1).'; 
tol = 1e-8; 
j = 1;
% plus self
r1 = 2*r; 
r2 = r+r; 
pass(j) = max(abs(r1(x)-2*fa(x))) < 1.5*tol;  % j = 1
j = j+1; 
pass(j) = max(abs(r2(x)-r1(x))) < 1e-15;  % j = 2
j = j+1; 

% plus constant: 
r1 = r+10; 
r2 = rfun(fa(x)+10, 'tol', 1e-9); 
pass(j) = max(abs(r2(x)-r1(x))) < tol;  
j = j+1; 

%minus: 
r1 = -r; 
pass(j) = max(abs(r1(x)+r(x)))==0;
j = j+1;

r1 = .5*r; 
r2 = r-r1; 
pass(j) = max(abs(r1(x, 'values')-.5*fa(x)))<5*tol;
j = j+1; 
r1 = r-r;
pass(j) = all(abs(r1(x))==0);

%plus two rfuns: 
[rb, fb] = gallery_rfun('wild');
h = @(x) fb(x) + fa(x); 
h2 = r + rb; 
pass(j) = max(abs(h2(x)-h(x))) < 5e1*tol;  

%plus an rfun and efun:


end