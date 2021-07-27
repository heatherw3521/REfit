function val = get( r, propName )
% get method for rfuns.
%   val = get(r, propName) returns the property specified in the string 
%   propName.Valid entries for the string propName are:
%    'nodes' = barycentric nodes/interpolating points/support points
%    'weights' = barycentric weights
%    'vals' = value of r at the nodes
%    'poles'  = poles of trigonometric barycentric rational
%    'residues' = residues assoc w the poles
%    'const' = additive constant assoc. with r.
%    'domain'  
%    'tol'   = tolerance parameter from original fit to data
%    'res'   = # of points in the sample from which r was constructed.
%%

% Get the properties.
switch ( propName )
    case 'domain'
        val = r.domain;
    case 'nodes'
        val = r.nodes;
    case 'vals'
        val = r.vals;
    case 'weights'
        val = r.weights;
    case 'poles'
        val = r.poles; 
    case 'residues'
        val = r.residues; 
    case 'const'
        val = r.const;
    case 'scl'
        val = r.scl;
    case 'tol'
        val = r.tol; 
    case 'res'
        val = r.res;
    otherwise
        error('RFUN:get:propName', ...
            [propName,' is not a valid property.'])
end

end
