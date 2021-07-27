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

k = 1; 
type = 'handle'; 
%parse input: 
for j = 1:nargin-1
    up = varargin{j};
    if isnumeric(up)
        k = up;
    else
        type = up; 
    end
end

%deal with constants and zero function: derivative is zero function.
if all(r.vals==0)
   R = rfun([0;0], [r.domain(1)+1/diff(r.domain);...
        r.domain(1)+(diff(r.domain)-1)/diff(r.domain)]);
   R.domain =r.domain;
   R.const = 0;
   R.scl = 1;
   h = @(x) 0*x;
   if strcmpi(type,'handle')
       varargout{1} = h;
       return
   elseif strcmpi(type, 'efun')
       R = ft(R); %converts to efun
   elseif strcmpi(type,'chebfun')
       R = chebfun(@(x) 0*x, r.domain); 
   end
   varargout{1} = R; 
   if nargout > 1
       varargout{2}= h;
   end
   return
end
       
%for handle, use the recursive formula from ():
h = @(x) diff_eval(r, x, k); 
if strcmpi(type, 'handle') %just return handle. 
    varargout{1} = h; 
    return
end

%deal with various return types: 
if strcmpi(type, 'chebfun')
    R = chebfun(h, r.domain, 'tol', r.tol, 'trig'); 
else %for rfuns and efuns: here it is easier to use efun than 
% to try to adaptively determine the appropriate grid.
    s = ft(r); 
    R = diff(s, k, 'efun'); %if type = efun, do nothing.
    if strcmpi(type, 'rfun')
        R = ift(R); %go from efun to rfun
    end
end
varargout{1} = R; 
if nargout >1
    varargout{2} = h;
end

end

%for now we use Fourier transforms and efun's diff. 
% this isn't the best practice since derivatives aren't always well
% represented with efuns. The closed-form formula will soon replace this:

function vals = diff_eval(r, x, k)
s = ft(r); 
ds = diff(s, k, 'efun'); 
vals = feval(ds, x, 'values');
end

% TO DO: code up closed-form formula and use this for function handle
% instead. The code is started below...

% function vals = diff_eval(r, x, k)
% %convert from [a, b] to [0, 1): 
% a = r.domain(1); 
% b = r.domain(2); 
% r.nodes = (r.nodes-a)./(b-a); 
% [a, b] = size(x); 
% x = x(:); 
% x = (x-a)./(b-a);
% r.nodes = zj; 
% r.weights = wj; 
% r.vals = fj; 
% % we compute d/dt r(x(t)), where t on [a, b] and dx/dt is 1/(b-a); 
% % the formula are recursive:
% %
% % pull out x that corresponds with zj and deal with them:
% [~, idxx, idxzj] = intersect(x, zj); 
% 
% if k==1
%         num = csc(pi*(x- zj.')).^2.*(f.'-feval(r, x)); 
%         num = num*wj;  
%         den = cot(pi*(xx- zj.'));
%         den = den*wj;
%         vals = -pi/(b-a)*num./den; 
%         if ~isempty(idxzj)
%             %compute vals at zj points:
%             zjj = zj(idxzj); 
%             vzj = zjj; %to do
%             vals(idxx) = vzj;
%         end
% end
    
