function [mins, maxs] = extrema(r)
% find locations for all local minima/maxima of r on its domain. 
% 
% 
%%
% see also rfun/min, rfun/max,rfun/roots.

% To find extrema, we differentiate, build an rfun, find its roots:
%rd.tol = min(1e-4, 1e2*r.tol); 
%

% TO DO: deal with strong singularities more appropriately.


rd = diff(r,1, 'rfun'); 
ext = roots(rd); ext = sort(ext); 
%h = diff(rd); %handle for second derivative (once we have 
% closed form formula in, this will be a very fast way to check.)

%for now, we can just check on a grid: 
L = length(ext);
mins = []; 
maxs = []; 
df = diff(ext); 
%extrema on endpoints: 
% LEFT
if abs(ext(1)-r.domain(1))<1e-3
    tr = ext(1) + .5*df(1);
    if (feval(r, tr)-feval(r, ext(1)))> 0
        mins = [mins; ext(1)]; 
        ext(1) = []; 
        df(1) = []; 
    elseif (feval(r, tr)-feval(r, ext(1)))< 0
        maxs = [maxs; ext(1)]; 
        ext(1) = []; 
        df(1) = []; 
    else %inflection point
        ext(1) = []; 
        df(1) = [];
    end
end

% RIGHT
if abs(ext(end)-r.domain(2))<1e-3
    tl = ext(1) - .5*df(end);
    if (feval(r, tl)-feval(r, ext(end)))> 0
        mins = [mins; ext(end)]; 
        ext(end) = [];  
    elseif (feval(r, tl)-feval(r, ext(1)))< 0
        maxs = [maxs; ext(end)]; 
        ext(end) = [];  
    else %inflection point
        ext(end) = []; 
    end
end
        
%check extrema: 
ex = [r.domain(1); ext; r.domain(end)]; 
df = diff(ex); 

for j = 2:length(ex)-1
    tl = ex(j)-df(j-1)/2; 
    tr = ex(j) + df(j)/2; 
    if (feval(r, tl)-feval(r, ex(j)))< 0 && (feval(r, tr)-feval(r, ex(j)))<0 
        maxs = [maxs; ex(j)]; 
    elseif (feval(r, tl)-feval(r, ex(j)))> 0 && (feval(r, tr)-feval(r, ex(j)))> 0
        mins = [mins; ex(j)];
    end 
end
      
end    

    



