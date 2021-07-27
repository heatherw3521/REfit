function out = residues(r)
% residues(s) returns the residues of the poles of r. 
%%

if isempty(r)
    out = []; 
end

out = r.residues;

end