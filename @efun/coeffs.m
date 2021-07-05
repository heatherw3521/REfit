function y = coeffs(s, varargin)
% get the Fourier coefficients of the trigonometric rational associated
% with s (given by ift(s) ). 
%
% getcoeffs(s, modes) returns the Fourier coefficients for the supplied
% modes. If modes are not specified, modes = [-s.res, s.res], where s.res
% is the bandlimit of the sample of r(x) used to construct s. 
%
% See also: efun\feval
%%

if isempty(s)
    y = []; 
    return
end

if isempty(varargin)
    a = s.res; 
    y = feval(s, (-a:a).'); 
else
    modes = varargin{1}; 
    y = feval(s, modes); 
end

end


