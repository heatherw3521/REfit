function pass = test_ift
warning off
%test the basic functionality of the ift function. 
j = 1; 
tol = 5e-6; 

% NOTE 1: the 'wild' function is one that efun struggles to resolve. 
% Using it to test the ift initiates the contigency loops in the ift, 
% which is what we want to test. (a better way to construct 'wild' is to 
% use rfun; then the ft if one wants an efun). 
%
% NOTE 2: a strange feature of the randomized Hankel approach + the fact
% that exponential sums are "sloppy models" is that
% the nodes/weights are not the same each time efun is called. 
% (maybe we should fix this by seeding ?). 

[s, fa] = gallery_efun('wild'); 
x = linspace(0, 1, 3000).'; 
%test polres: 
r = ift(s, 'polres'); 
pass(j) = max(abs(r(x) - fa(x)))< tol;
j = j+1; 
%test z-transform polres: 
z = exp(2*pi*1i*x); 
r2 = ift(s, 'zt'); 
pass(j) = max(abs(r2(z) - fa(x)))< tol;
j = j+1; 

%test efun to rfun via ift: 
r2 = ift(s); 
pass(j) = max(abs(r2(x) - fa(x)))< tol;
j = j+1; 
pass(j) = (length(r2)==2*length(s)); %correct # of poles?
j = j+1; 

%test aaa-version of ift: 
r2 = ift(s, 'aaa');
pass(j) = max(abs(r2(x) - fa(x)))< 1e1*tol;
j = j+1; 

%test inputting various parameters: 
r2 = ift(s, 2000); 
pass(j) = max(abs(r2(x) - fa(x)))< tol;
j = j+1; 

v = linspace(0, 1, 3000); 
r2 = ift(s, v); 
pass(j) = max(abs(r2(x) - fa(x)))< tol;
j = j+1;

r2 = ift(s, 'aaa', 2000); 
pass(j) = max(abs(r2(x) - fa(x)))< 1e3*tol;
j = j+1; 

r2 = ift(s, 'fft'); 
pass(j) = max(abs(r2(x) - fa(x)))< 1e1*tol;
j = j+1; 

r2 = ift(s, 'ifft', 2000); 
pass(j) = max(abs(r2(x) - fa(x)))< 1e1*tol;
j = j+1; 

r2 = ift(s, 'aaa', 'ifft'); 
pass(j) = max(abs(r2(x) - fa(x)))< 1e3*tol;
j = j+1; 

% a few tests with nonstandard domain:
% why these perform more poorly is something to further
% investigate. they should be the same
xx = linspace(3, 5, 3000).';  
x1 = linspace(0, 1, 2*4000+2).'; x1 = x1(1:end-1);
ss = efun(fa(x1), 'tol', 1e-7, 'domain', [3, 5]);  
r = ift(ss, 'polres');
pass(j) = max(abs(r(xx) - fa(x)))< 1e1*tol;
j = j+1; 

r2 = ift(ss); 
pass(j) = max(abs(r2(xx) - fa(x)))< 1e3*tol;
j = j+1; 
pass(j) = all(real(r2.poles) >= 3 & real(r2.poles) <= 5);
j = j+1; 
 
r = ift(ss, 'aaa');
pass(j) = max(abs(r(xx) - fa(x)))< 1e3*tol;
j = j+1; 
pass(j) = all(real(r2.poles) >= 3 & real(r2.poles) <= 5);
 
end

