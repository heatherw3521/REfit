function varargout = sample(s, varargin)
% sample an efun or its inverse Fourier transform. 
%
% [vals, loc] = sample(s) returns a sample (vals) of r(x) = ift(s) on
% its domain on the sampling grid (loc) that was used to construct s.
%
% [vals, loc] = sample(s, 'Fourier'). Samples s in Fourier space. 

if isempty(s)
    if nargout==2
        varargout = {[],[]}; 
    else
        varargout = {[]};
    end
    return
end

res = s.res; 
Fourier_flag = 0;
for j = 1:nargin-1
    up = varargin{j}; 
    if isnumeric(up)
        res = up; 
    else %up = 'Fourier'. Get coeffs.
        Fourier_flag = 1;
    end
end
dom = s.domain; 
if ~Fourier_flag
    %sample in time: 
    loc = linspace(dom(1), dom(2), 2*res+2); loc = loc(1:end-1).';
    vals = feval(s, loc, 'values'); 
else
    %sample in Fourier space: 
    loc = (-res:res).';
    vals = feval(s, loc); 
end

if nargout==2
    varargout = {vals, loc}; 
else
    varargout = {vals};
end

end
    
    