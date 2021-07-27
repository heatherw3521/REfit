function varargout = cumsum(s, varargin)
% indefinite integration of an efun. 
%
% cumsum(s) returns a function handle for computing 
% the Fourier coefficients of the function 
% F(x) = integral( ift(s)-mean(ift(s)) over [a, x] 
% where ift = inverse Fourier transform, and [a, b] 
% is the domain of ift(s). 
% 
% cumsum(s, 'type') returns F(x) as an object based on the input 'type'. 
% One can choose 'type' = 'rfun', 'efun', or 'chebfun'.  
%
%%
% See also, efun\sum, efun\integral

%%
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

%% 
type = 'handle'; 
%%
% if ift(s) is not a mean zero function, we subtract off the mean.
if abs(s.const)>1e-13
    warning('efun:cumsum: cumsum is performed on efuns representing mean zero functions. Setting mean to zero now.') 
end
s.const = 0; 

%parse input
if ~isempty(varargin)
    type = varargin{1}; 
end
        
%%
% build function handle: 
h = @(j) eval_coeff(s, j);
d = s.res; 

if strcmpi(type, 'handle') %wants function handle
    varargout = {h};
    return
else                       %wants efun, rfun, or chebfun
    S = efun( h((0:d)'), 0:d, 'coeffs', 'tol',1e-10); 
end
S.domain = s.domain;

if strcmpi(type, 'rfun')
    S = ift(S);
elseif strcmpi(type, 'chebfun')
    S = chebfun(@(x) feval(S, x, 'values'), s.domain, 'tol', S.tol); 
end
 
if nargout > 1
    varargout = {S, h};
else
    varargout = {S}; 
end

%%
end

function vals = eval_coeff(s, j)
%s is efun, j = coeffs to eval
dom = s.domain; 
L = dom(2)-dom(1); 
j = j(:);  
vals =  L*feval(s, j)./(2*pi*1i*j); 
vals((j==0)) = 0; %fix where j = 0; 
end

