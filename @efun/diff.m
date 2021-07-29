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
%%
% See also rfun/diff

if isempty(s)
    varargout{1} = []; 
    if nargout > 1
        varargout{2} = [];
    end
    return
end

if strcmp(s.space, 'time')
    error('efun:diff: differentiation for efuns in signal space not available.')
end

%set default:
type = 'handle'; 
k = 1; 

%parse input:
if ~isempty(varargin)
    len = nargin-1;
    for jj = 1:len
        up = varargin{jj}; 
        if isnumeric(up)
            k = up; 
        else
            type = up; 
        end
    end
end
       
h = @(j) eval_coeff(s, j, k); 
d = s.res; 

if strcmpi(type, 'handle') %wants function handle
    varargout = {h};
    return
else                       %wants efun, rfun, or chebfun
    M = min(8000, 10*d); 
    S = efun( h((0:M)'), 0:M, 'coeffs', 'tol',1e-10); 
end
S.domain = s.domain;

if strcmpi(type, 'rfun')
    S = ift(S);
elseif strcmpi(type, 'chebfun')
    S = chebfun(@(x) feval(S, x, 'values'), s.domain, 'tol', S.tol, 'trig'); 
end
 

if nargout > 1
    varargout = {S, h};
else
    varargout = {S}; 
end
    
end

function vals = eval_coeff(s, j, k)
j = j(:); 
P = s.domain(2)-s.domain(1); 
vals =  ((2*pi*1i*j/P).^k).*feval(s, j); 
end
