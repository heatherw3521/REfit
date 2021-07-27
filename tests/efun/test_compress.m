function pass = test_compress()

j = 1; 
pass(j) = isempty(compress(efun())); 
j = j + 1; 

[s, ~] = gallery_efun('spline');
ss = compress(s, 1e-4); 

pass(j) = (length(ss) <= length(s));
j = j + 1; 

x = linspace(0, 1, 1000).'; 
pass(j) = max(abs(ss(x, 'values') - s(x, 'values'))) < 1e-3; 

end
