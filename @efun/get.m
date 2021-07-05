function val = get( s, propName )
%GET       GET method for efuns.
%   P = GET(s, PROP) returns the property P specified in the string PROP from
%   Valid entries for the string PROP are:
%    'exp'
%    'weights'
%    'domain'
%    'sv' % singular values of the Hankel matrix from Prony's method.
%    'space' % exp sum is in time domain or Fourier space?
%    'tol'   %tolerance parameter from original fit to data
%    'scl'   %scaling parameter
%    'const' %constant parameter
%    'res'   %bandlimit of signal that sum was originally fit to. 
%%

% Get the properties.
switch ( propName )
    case 'domain'
        val = s.domain;
    case 'exp'
        val = s.exp;
    case 'weights'
        val = s.weights;
    case 'sv'
        val = s.sv; 
    case 'space'
        val = s.space;
    case 'tol'
        val = s.tol;
    case 'scl'
        val = s.scl;
    case 'const'
        val = s.const;   
    case 'res'
        val = s.res;
    otherwise
        error('EFUN:get:propName', ...
            [propName,' is not a valid property.'])
end
end
