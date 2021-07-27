function pass = test_times()
%test the basic functionality of the times operator: 
tol = 1e-9; 
j = 1; 
[s, fa] = gallery_efun('spline'); 
r = ift(s); 
[f, x] = sample(r); 

r2 = 3*r; 
pass(j) = max(abs(r2(x)- 3*f))< tol; 
j = j+1; 

r2= 3.*r; 
pass(j) = max(abs(r2(x)- 3*f))< tol; 
j = j+1;

r2 = r.*r;
pass(j) = max(abs(r2(x)- f.*f))< tol; 
j = j+1; 

r1 = r; 
r1.domain = [-1, 1]; 
pass(j) = 0;
try
    R = r.*r1;
catch 
    pass(j) = 1;
end
j = j+1; 

%test times a zero
r1 = ift(efun([0;0])); 
r2 = r.*r1; 
pass(j) = all(r2(x)==0);
j = j+1; 

%test times a constant
r1.const = 5; 
r2 = r1.*r; 
pass(j) = max(abs(r2(x)-5*f))<tol; 
j = j+1;

%test rfun times function_handle
r2 = r.*fa; 
pass(j) = max(abs(r2(x)- f.*f))< 5*tol; 
j = j+1; 

%test efun times rfun
s = ft(r); 
r2 = r.*s; 
pass(j) = max(abs(r2(x)- f.*f))< tol;

end

