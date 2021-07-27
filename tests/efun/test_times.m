function pass = test_times()
%test the basic functionality of the times operator: 
tol = 1e-9; 
j = 1; 
[s, fa] = gallery_efun('spline'); 
[f, x] = sample(s); 

s2 = 3*s; 
pass(j) = max(abs(s2(x, 'values')- 3*f))< tol; 
j = j+1; 

s2= 3.*s; 
pass(j) = max(abs(s2(x, 'values')- 3*f))< tol; 
j = j+1;

s2 = s.*s;
pass(j) = max(abs(s2(0:100)- s(0:100).*s(0:100)))< tol; 
j = j+1; 

s1 = s; 
s1.domain = [-1, 1]; 
pass(j) = 0;
try
    S = s.*s1;
catch 
    pass(j) = 1;
end
j = j+1; 

%test efun times function_handle
f = @(j) s(j); 
s2 = s.*f; 
pass(j) = max(abs(s2(0:100)- s(0:100).*s(0:100)))< tol; 
j = j+1; 

%test efun times rfun
r = ift(s); 
s2 = s.*r; 
pass(j) = max(abs(s2(0:100)- s(0:100).*s(0:100)))< tol;

end

