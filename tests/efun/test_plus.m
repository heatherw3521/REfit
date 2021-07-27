function pass = test_plus()
%test plus function 
[s, fa] = gallery_efun('chirp');
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

%minus: 
s1 = -s; 
pass(j) = max(abs(s1(x, 'values')+s(x, 'values')))<tol;
j = j+1;

s1 = .5*s; 
s2 = s-s1; 
pass(j) = max(abs(s1(x, 'values')-.5*fa(x)))<5*tol;
j = j+1; 
s1 = s-s;
pass(j) = all(abs(s1(0:100))==0);



%plus two efuns: 
[g, fb] = gallery_efun('wild');
%h = efun(fb(x) + fa(x), 'tol', 1e-10);
h = @(x) fb(x) + fa(x); 
h2 = s + g; 
pass(j) = max(abs(h2(x, 'values')-h(x))) < 5e2*tol;  

%plus an rfun and efun:


end