function y = coeffs(r, varargin)
% get the Fourier coefficients of an rfun
%
% getcoeffs(r, modes) returns the Fourier coefficients for the supplied
% modes. If modes are not specified, modes = [-b, b], where b is determined
% based on (1) the size of the original sample r was fitted to and (2) the 
% decay of the Fourier coefficients.
%
% See also: rfun/ft, efun/coeffs, efun/feval. 
%%
if isempty(r)
    y = []; 
    return
end

s = ft(r); 
y = coeffs(s, varargin{:}); 

end


