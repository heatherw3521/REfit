function [res, const] = residues(s, varargin)
% residues(s) returns the residues of the poles of r(z), 
% where r = ift(s), r(z) = r(x) on the unit circle,
% and z = exp(2*pi*1i*(x-a)/(b-a)) with x defined on [a, b].


 if isempty(s)
    res = [];
    const = []; 
    return
 end

 if strcmpi(varargin{1}, 'zt')
    w = s.weights; 
    const = s.const; 
    scl = s.scl; 
    res = scl*[w.*exp(s.exp); -conj(w).*exp(-conj(s.exp))];
 else
    error('for now residues are only available with zt flag')
 end

end