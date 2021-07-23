function y = coeffs(r, varargin)
% get the Fourier coefficients of an rfun
%
% getcoeffs(r, modes) returns the Fourier coefficients for the supplied
% modes. If modes are not specified, modes = [-b, b], where b is 
% adaptively determined so that for abs(m) > b, abs(m) < eps. 
%
% See also: rfun/ft, efun/coeffs, efun/feval. 
%%
if isempty(r)
    y = []; 
    return
end

s = ft(r); % this will automatically determine approx. bandlimit.
y = coeffs(s, varargin{:}); 

end


