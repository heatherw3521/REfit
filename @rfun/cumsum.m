function varargout = cumsum(r, varargin)
% indefinite integration of an rfun. 
%
% cumsum(r) returns a function handle for computing 
% the integral of r over [r.domain(1), x].  
% 
% cumsum(r, 'type') returns F(x) as an object based on the input 'type'. 
% One can choose 'type' = 'rfun', 'efun', or 'chebfun'.  For efun and rfun
% types, F(x) must be periodic, so it is assumed that mean(r)=0. 
%
% See also, rfun\sum, rfun\integral

%%
if isempty(r)
    varargout{1} = []; 
    if nargout > 1
        varargout{2} = []; 
    end
    return
end

a = r.domain(1); 

type = 'handle';
if ~isempty(varargin)
    type = varargin{1}; 
end
%for now, call efun. Later we will add option to do GL quadrature.

%handle the constant term: 
const = r.const; 
r.const = 0; 
s = cumsum(ft(r), 'efun'); 
d = r.res; 
h = @(x) eval_cumsum(s, x, const, a);
if strcmpi(type, 'handle') %wants function handle
    varargout = {h};
    return
elseif strcmpi(type, 'rfun')
    if abs(const) >1e-12 %not a mean zero function
        warning('rfun:cumsum: setting mean to zero.')
    end
    S = ift(s);
elseif strcmpi(type, 'efun')
    if abs(const) >1e-12 %not a mean zero function
        warning('rfun:cumsum: setting mean to zero.')
    end
    S = s; 
elseif strcmpi(type, 'chebfun')
    S = chebfun(h, r.domain, 'tol', r.tol); 
end
 
if nargout > 1
    varargout = {S, h};
else
    varargout = {S}; 
end


%%

end

function val = eval_cumsum(s, x, const,a)
    val = feval(s,x, 'values') + (x-a)*const;
end



