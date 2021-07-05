function [coeffs, fvals, s_grid] = get_sample(f, space, dom, tol) 
% adaptively selects a sample f on domain dom. 
%
% the size of the sample depends on tol and uses a variation 
% of Chebfun's "chop" function. 
%
% see: (LITERATURE)

happy = 0; 
s_grid = (dom(2)-dom(1))*(0:4).'/5 + dom(1); %initialize grid.
fvals = f(s_grid); 

while ~happy
    %double the grid and sample new points. 
    dist = (s_grid(2)-s_grid(1))/2;
    s_grid_old = s_grid; 
    s_grid = zeros(2*length(s_grid_old),1); 
    s_grid(2:2:end) = s_grid_old + dist; 
    s_grid(1:2:end) = s_grid_old; 
    fvals_old = fvals; 
    fvals = zeros(2*length(fvals_old),1); 
    fvals(1:2:end) = fvals_old; 
    fvals(2:2:end) = f(s_grid(2:2:end)); 
    %get coeffs and check tail
    coeffs = sample2coeffs(fvals, 'pos'); 
    cutoff = standardChop(coeffs, tol); 
    %are we happy yet?
    happy = ~(cutoff == length(coeffs));
    
    if length(coeffs) >= 20000
        error('efun:constructor:get_sample:sufficient sampling is not possible.')
    end    
end

if strcmpi(space, 'time')
    coeffs = fvals; 
end

end