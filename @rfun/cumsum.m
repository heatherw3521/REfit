function varargout = cumsum(s, varargin)
% indefinite integration of an efun. 
%
% cumsum(s,c) returns a function handle for computing 
% the integral of r over [c, x] \subset [a b], where [a, b] 
% is the domain of ift(s). By default, c = a. 
% 
% cumsum(s, 'type') returns F(x) as an object based on the input 'type'. 
% One can choose 'type' = 'rfun', 'efun' (returns ift(cumsum(r)), or 'chebfun'.  
%
% See also, rfun\sum, rfun\integral

%%
% for now, we only consider efuns representing Fourier series
if strcmp(s.space, 'time')
    error('efun:diff: differentiation for efuns in signal space not available.')
end

%%
[a, ~] = r.domain; 

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

% what is the best way to do this? (1) get an efun and use ft/nufft to
% evaluate? (2) GL quadrature?

end


