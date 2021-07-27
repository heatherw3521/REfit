function out = residues(s, varargin)
% residues(s, 'zt') returns the residues of the poles of r(z), 
% where r = ift(s), r(z) = r(x) on the unit circle,
% and z = exp(2*pi*1i*(x-a)/(b-a)) with x defined on [a, b].

if isempty(s)
    out = []; 
end

if isempty(varargin)
    out = 'to do';
else
    out = w.*exp(s.nodes); 
end

end