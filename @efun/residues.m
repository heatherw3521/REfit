function out = residues(s, varargin)
% residues(s) returns the residues of the poles of r(x) = ift(s). 
%
% residues(s, 'zt') returns the residues of the poles of r(z), 
% where r(z) = r(x) on the unit circle and z = exp(2*pi*1i*x). 

if isempty(s)
    out = []; 
end

if strcmp(s.space, 'value')
    error('efun:getpoles:cannot get residues for an efun in time/signal space')
end

if isempty(varargin)
    out = 'to do';
else
    out = w.*exp(s.nodes); 
end

end