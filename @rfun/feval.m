function vals = feval(r, zz) 
% Evaluate an rfun R at the points zz. 
%
% See also: rfun/sample, rfun/coeffs.

[s1, s2] = size(zz); 
zz = zz(:);
dom = r.domain; 
a = dom(1); 
b = dom(2); 
%map to [0, 1): 
zz = (zz-a)/(b-a); 
zj = (r.nodes-a)/(b-a); 

fj = r.vals; 
wj = r.weights; 
const = r.const; 
scl = r.scl; 

% Evaluate rational function in barycentric form.

CC = bsxfun(@(vals, params) cot(pi*(vals - params)), zz, zj.'); 
vals = (CC*(wj.*fj))./(CC*wj);             % vector of values

% Deal with infs: r(inf) = lim r(zz) = 0, since r is type (n-1, n)
%vals(isinf(zz)) = sum(wj.*fj)./sum(wj);
vals(isinf(zz)) = 0; 

% Deal with NaN:
ii = find(isnan(vals));
for jj = 1:length(ii)
    if ( isnan(zz(ii(jj))) || ~any(zz(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at interpolating points.
        vals(ii(jj)) = fj(zz(ii(jj)) == zj);
    end
end

% Reshape to input format:
vals = scl*vals+const; 
vals = reshape(vals, s1,s2);

end 
    