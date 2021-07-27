function R = compress(r, varargin)
% compress an rfun. This function tries to represent r with a 
% trigonometric rational function of smaller degree type. 
%
% tolerance can be specified with compress(r, tol). 
% Default tol = r.tol. 
%
%%
% See also: rfun/ft, efun/compress.

%%
%%empty case:
if isempty(r) 
    %pass to constructor 
    R = rfun(); 
    return
end

m = length(r);
%length 2 case: no compression is possible
if (m==2) 
    R = r; 
    return
end

tol = r.tol; 
if ~isempty(varargin)   
    tol = varargin{1};   
end

%%
% call efun via Fourier transform to compress: 
s = ft(r, 'tol', tol); 
% construct an rfun: 
R = ift(s); 

end
    






