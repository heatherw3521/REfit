function varargout = diff(r, varargin)
% differentiation of an rfun 
%
% diff(r, k) returns a function handle for computing the  
% the kth derivative of r. By default, k = 1.
%
% diff(s, k, 'type') returns the kth derivative of r as an
% object based on the input 'type'. One can choose
% 'type' = 'rfun', 'efun'(returns ift(diff(r))), or 'chebfun'. 
%
% See also efun/diff

if isempty(r)
    R = []; 
    return
end

len = length(varargin);
if len == 1
    if isnumeric(varargin{1})
        k = varargin{1}; 
    else
        type = varargin{1};
        k = 1; 
    end
elseif len ==2
    k = 1; 
    type = varargin{2}; 
else
    k = 1; 
    type = 'handle'; 
end


end

function vals = eval_diff(r, j, k)

end
