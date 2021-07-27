function varargout = coeffs2vals(varargin)
%wrapper for coeffs2sample
if nargout==2
    [varargout{1}, varargout{2}] = coeffs2sample(varargin{:});
else
    varargout{1} = coeffs2sample(varargin{:}); 
end
