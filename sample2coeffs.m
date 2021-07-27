function F = sample2coeffs(f, varargin)
% gets correct vector of Fourier coeffs to pass to efun constructor.
%
% f can be a vector or a function handle or chebfun. 
% With function handle:
% supply domain = [a, b] and number of coeffs reqd. 
%
% with chebfun: supply # coeffs (sample rate)
% 
% if f is a vector, it should be a vector of equispaced values
% of f sampled over domain. 

if isa(f, 'function_handle')
    dom = varargin{1}; 
    N = varargin{2}; 
    x = linspace(dom(1), dom(2), N+1); 
    x = x(1:end-1); 
    f = f(x); 
elseif isa(f, 'chebfun')
    dom = [f.domain(1), f.domain(end)];
    N = varargin{1};
    x = linspace(dom(1), dom(2), N+1); 
    x = x(1:end-1); 
    f = f(x); 
end
f = f(:); 
f = 1/length(f)*fft(f);
F = fftshift(f); 

if (any(strcmpi(varargin, 'pos')) || any(strcmpi(varargin, 'positive'))...
|| any(strcmpi(varargin, 'nonneg')))
 
%get the non-neg indexed coeffs only:
    n = length(F); 
    F = F((n-mod(n,2))/2+1:end);
end
