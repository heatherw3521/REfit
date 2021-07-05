function varargout = subsref(s, index)
%subsref for efuns. 
%( )
%s(x), where x is a matrix, returns the 
% matrix of associated values of s at x. 


idx = index(1).subs;
switch index(1).type
    
case '()' %eval
    x = idx{1};
    
   [a, b] = size(x); 
    %flatten for eval: 
    x = x(:);    
    out = feval(s, x, idx{2:end});
    if ~(a==1) && ~(b ==1) %matrix needs reshaped
        out = reshape(out, a, b);  
    end
    varargout = { out };
    
 case '.' %get property
     out = get(s, idx);
     varargout = { out };   
     
     
end
 

end
