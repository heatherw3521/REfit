function varargout = cumsum(s, varargin)
% indefinite integration of an efun. 
%
% cumsum(s,c) returns a function handle for computing 
% the Fourier coefficients of the function 
% F(x) = integral( ift(s)-mean(ift(s)) over [c, x] \subset [a b], 
% where ift = inverse Fourier transform, and [a, b] 
% is the domain of ift(s). By default, c = a. 
% 
% cumsum(s, 'type') returns F(x) as an object based on the input 'type'. 
% One can choose 'type' = 'rfun', 'efun', or 'chebfun'.  
%
% See also, efun\sum, efun\integral

%%
% for now, we only consider efuns representing Fourier series
if strcmp(s.space, 'time')
    error('efun:diff: differentiation for efuns in signal space not available.')
end

%%
[a, ~]= s.domain; 
%%
% if ift(s) is not a mean zero function, we subtract off the mean. 
s.const = 0; 

%parse input
if isempty(varargin)
    type = 'handle';
    c = a;
elseif length(varargin)==1
    if isa(varargin{1}, 'string')
        type = varargin{1};
        c = a; 
    else
        c = varargin{1};
        type = 'handle';
    end
elseif isa(varargin{1}, 'string')
    type = varargin{1}; 
    c = varargin{2}; 
else
    type = varargin{2}; 
    c = varargin{1}; 
end

%%
% build function handle: 
h = @(j) (s.scl)*eval_coeff(s, j, k);
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

