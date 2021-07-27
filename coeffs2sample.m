function varargout = coeffs2sample(coeffs, varargin)
% transforms a sample of coefficients over modes [-M:M] (or
% modes [0:M] with flag 'pos') to values in time on an 
% equispaced grid. 
%
% one can include 'domain', [a, b] to set the domain to [a, b]. 
% By default, [a, b] = [0, 1]. 
% 
% See also sample2coeffs.

dom = [0, 1]; 
pos_flag = false; 
if ~isempty(varargin)
    dom_id= find(strcmpi(varargin, 'domain'));
    if ~isempty(dom_id)
        dom = varargin{dom_id+1}; 
    end
    if any(strcmpi(varargin, 'pos'))
        pos_flag= true; 
    end
end

if pos_flag
    n = length(coeffs); 
    coeffs = [flip(conj(coeffs(2:end)));coeffs]; 
else
    n = (length(coeffs)-1)/2; 
end
coeffs = ifftshift(coeffs); 
vals = length(coeffs)*ifft(coeffs); 
vals = real(vals); 
pts = linspace(dom(1),dom(2), 2*n); 
pts = pts(1:end-1).';

varargout{1} = vals; 
if nargout==2
    varargout{2} = pts; 
else
    varargout{1} = vals; 
end
end