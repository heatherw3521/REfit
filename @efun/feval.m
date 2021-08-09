function vals = feval(F, x, varargin) 
% evaluate an efun F at the points in x. 
% 
% feval(F, x) evaluates the exponential sum at locations x. 
%
% feval(F, x, 'time')or feval(F, x, 'values')  evaluates the 
% inverse Fourier transform of F at locations x. 
%
% feval(F, z, 'zt') evaluates the rational r associated with F in the 
% complex z-plane. F is the Fourier transform of r(x), where r(x) = r(z) 
% for z restricted to the unit circle and z = exp(2*pi*1i*(x-a)/(b-a)), 
%  x on [a, b]. 
%
%%
% See also efun/plot, efun/ift. 
 
sz = size(x); 
x = x(:); 
r = F.exp; 
w = F.weights; 
scl = F.scl; 
const = F.const; 

%if pr, rat, value keywords are present, evaluate in PR form: 
keys = {'pr', 'polres', 'rat', 'rational', 'value', 'values', 'time'}; 
if strcmpi(F.space, 'Fourier') && ~isempty(varargin)
    evaltype = varargin{1};
    if any(strcmpi(evaltype, keys))
        h = ift(F, 'polres'); 
        vals = h(x); 
        return
    end
    if any(strcmpi(evaltype, 'zt'))
        h = ift(F, 'zt'); 
        vals = h(x); 
        return
    end
end
%evaluate exponential sum directly: 

%separate out neg and non-neg:
if strcmpi(F.space, 'Fourier')
    [idx1, ~] = find(x < 0);
    x1 = x(idx1); 
    [idx2, ~] = find(x >= 0); 
    x2 = x(idx2); 
else %for exponential sums in value space, we directly eval all points
    x2 = x; 
    x1 = []; 
    idx2 = 1:length(x); 
    idx1 = []; 
end

vals2 = []; 
if ~isempty(x2)
    M = repmat(r.', length(x2),1); 
    M = M.*(repmat( x2, 1,length(r)));
    M = exp(M);  
    vals2 = scl*(M*w); 
    if strcmpi(F.space, 'Fourier')
        vals2(x2==0) = const; %fix the zero-coeff value. 
    end
end

vals1 = []; 
if ~isempty(x1)
    x1 = abs(x1); 
    M = repmat(r.', length(x1),1); 
    M = M.*(repmat( x1, 1,length(r)));
    M = exp(M);  
    vals1 = scl*(M*w);
    vals1 = conj(vals1); 
end
vals = zeros(length(x),1); 
vals(idx1)= vals1; 
vals(idx2) = vals2; 

%if in time/value/signal space, the values should always be real: 
if strcmpi(F.space, 'time')
    vals = real(vals + const);     
end

vals = reshape(vals, sz); 
end