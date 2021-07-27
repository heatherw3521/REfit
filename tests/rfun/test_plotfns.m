function pass = test_plotfns
%test various plot functions for efuns: 

r = rfun(@(x) exp(sin(2*pi*(x+1)))); 
j = 1; 

%plot: 
try
    plot(r);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1; 

r.domain = [-1, 1]; 
try
    plot(r);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

%plotcoeffs: 
r.domain = [0,1];
try
    plotcoeffs(r);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

%poleplot: 
try
    poleplot(r);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

try
    poleplot(r, 'zt');
    pass(j) = true;
catch ME 
    pass(j) = false;
end

end




