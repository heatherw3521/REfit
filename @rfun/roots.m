function zer = roots(r, varargin)
% Find the roots of r. 
% 
% roots(r) finds the real-valued roots of r(x) on its domain. 
% roots(r, 'all') finds all roots, including complex valued ones. 
%
% See also: rfun/roots

%% Compute poles, residues and zeros of the rational.

wj = r.weights; 
zj = r.nodes; 
fj = r.scl*r.vals; %interpolation condition with constant 
dom = r.domain; 

m = length(wj);
a = dom(1); 
b = dom(2); 
% Compute poles via generalized eigenvalue problem:
B = diag([ones(m, 1); 0]); 
B(1:m, end) = -1i*wj; 
E = diag([exp(1i*2*pi/(b - a)*zj);0]); 
E(1:m,end) = 1i*wj.*exp(1i*2*pi/(b - a)*zj); 
E(end, 1:m) = fj+r.const; 

% add refinement for real roots, what about complex roots? 


zer = eig(E, B);
%remove zeros at z = \pm infinity and z = 0:
zer = zer(abs(zer)> 1e-6); 
zer = zer((abs(zer)<=1e10)); 
zer = zer(~isinf(zer));
%put in terms of x
zer = log(zer)/(2*pi/(b - a)*1i);

%correct zeros so that they are on the interval [a, b]: 
zer = a + mod(real(zer), b-a) + 1i*imag(zer); 

if isempty(varargin)
    zer = zer(abs(imag(zer))<1e-13); 
    zer = real(zer); 
end

    
end
