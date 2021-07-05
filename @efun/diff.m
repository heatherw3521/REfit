function varargout = diff(s, varargin)
% differentiation of an efun 
%
% diff(s, k) returns a function handle for computing the  
% the Fourier coeffs of the kth derivative of r(x) = ift(s). 
% By default, k = 1.
%
% diff(s, k, 'type') returns the kth derivative of ift(s) as an
% object based on the input 'type'. One can choose
% 'type' = 'rfun', 'efun', or 'chebfun'. 
%
% See also rfun/diff

if isempty(s)
    S = []; 
    return
end

if strcmp(s.space, 'time')
    error('efun:diff: differentiation for efuns in signal space not available.')
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

h = @(j) (s.scl)*eval_coeff(s, j, k); 
d = s.res; 

if strcmpi(type, 'handle') %wants function handle
    varargout = {h};
    return
else                       %wants efun, rfun, or chebfun
    S = efun( h((0:10*d)'), 0:10*d, 'coeffs', 'tol',1e-10); 
end
S.domain = s.domain;

if strcmpi(type, 'rfun')
    S = ift(S);
elseif strcmpi(type, 'chebfun')
    S = chebfun(@(x) feval(S, x, 'values'), s.domain, 'tol', S.tol, 'trig'); 
end
 

if length(varargout) > 1
    varargout = {S, h};
else
    varargout = {S}; 
end
    
 

end

function vals = eval_coeff(s, j, k)
j = j(:); 
vals =  (2*pi*1i*j).^k.*feval(s, j); 
end
