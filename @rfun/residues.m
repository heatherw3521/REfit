function out = residues(r, varargin)
% residues(s) returns the residues of the poles of r. 
%
% residues(s, 'zt') returns the residues of the poles of r(z), 
% where  z = exp(2*pi*1i*x/(b-a)). 

if isempty(r)
    out = []; 
end


if isempty(varargin)
    out = r.residues;
else
    out = 'to do'; 
end

end