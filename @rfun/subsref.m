function varargout = subsref(r, index)
%subsref for rfuns. 
%( )
%r(x), where x is a vector, returns the 
% vector of associated values of r at x. 

idx = index(1).subs;
switch index(1).type
    
case '()' %eval
    x = idx{1};
    
   [a, b] = size(x); 
    %flatten for eval: 
    x = x(:);    
    out = feval(r, x);
    if ~(a==1) && ~(b ==1) %matrix needs reshaped
        out = reshape(out, a, b);  
    end
    varargout = { out };
    
 case '.' %get property
     out = get(r, idx);
     varargout = { out };   
 end
end
