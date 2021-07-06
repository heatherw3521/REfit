function pass = test_plotfns
%test various plot functions for efuns: 

s = efun(@(x) exp(sin(2*pi*(x+1)))); 
j = 1; 

%plot: 
try
    plot(s);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1; 

s.domain = [-1, 1]; 
try
    plot(s);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

%plotcoeffs: 
try
    plotcoeffs(s);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

%poleplot: 
try
    poleplot(s);
    pass(j) = true;
catch ME 
    pass(j) = false;
end
j = j+1;

try
    poleplot(s, 'zt');
    pass(j) = true;
catch ME 
    pass(j) = false;
end

end




