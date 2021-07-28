function pass = test_ft()
%test the Fourier transform: 
tol = 1e-9; 
j = 1; 
x = linspace(0, 1, 4000).'; x = x(1:end-1);
[s, ~] = gallery_efun('spline'); 
r = ift(s); 
s2 = ft(r); 
pass(j) = length(r)>=2*length(s2); %1
j = j+1;
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < tol; %2
j = j+1; 


s2 = ft(r, 2000); %input a sample size parameter
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < tol; %4
j = j+1; 

s2 = ft(r, 'tol', 1e-5); %try passing tol parameter
pass(j) = max(abs(s(x, 'values')-s2(x, 'values'))) < 1e-5; %5
j = j+1;
r2 = ift(s2); 
pass(j) = max(abs(r2(x)-r(x))) < 1e-5; %6
j = j+1;

%try with different domain: 
s.domain = [3, 5]; 
x2 = linspace(3, 5, 4000).'; x = x2(1:end-1);
r = ift(s); 
s2 = ft(r); 
pass(j) = (length(r)>=2*length(s2)); %7
j = j+1;
pass(j) = max(abs(s(x2, 'values')-s2(x2, 'values'))) < 5*tol; %8
j = j+1; 

%try a different function:  
[r, ~]= gallery_rfun('wild');
s = ft(r); 
pass(j) = length(r)>=2*length(s); %9
j = j+1;
pass(j) = max(abs(s(x, 'values') - r(x))) < tol;  %10

end






