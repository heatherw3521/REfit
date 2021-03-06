function pass = test_sample
%test basic functionality of the sample function

[s, ~] = gallery_efun('spline'); 
n = s.res; 
x = linspace(0, 1, 2*n+2); 
x = x(1:end-1).';

%two outputs: 
[f, xx] = sample(s); 
j = 1; 
pass(j) = all(abs(xx - x)==0);
j = j + 1; 
pass(j) = all(abs(f - s(x, 'values'))==0);
j = j+1; 

%one output: 
ff = sample(s); 
pass(j) = all(abs(ff-f)==0); 
j = j+1; 

% Fourier space: 
[F, N] = sample(s, 'Fourier'); 
pass(j) = all(abs((-n:n).' - N)==0);
j = j+1; 
pass(j) = all(abs(F-s((-n:n).'))==0); 
j = j+1; 
%one output: 
F = sample(s, 'Fourier'); 
pass(j) = all(abs(F-s((-n:n).'))==0); 
j = j+1; 

%test empty: 
[a, b] = sample(efun()); 
pass(j) = all(isempty( [a,b]));
j = j+1; 
a = sample(efun); 
pass(j) =  isempty(a);


end







