function pass = test_compress()

j = 1; 
pass(j) = isempty(compress(rfun())); 
j = j + 1; 

[s, fa] = gallery_efun('spline');
r = ift(s); 
rr = compress(r, 1e-4); 

pass(j) = (length(rr) <= length(r));
j = j + 1; 

x = linspace(0, 1, 1000).'; 
pass(j) = max(abs(rr(x) - fa(x))) < 1e-3; 

end