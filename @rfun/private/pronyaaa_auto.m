function [r, pol, res, zj, fj, wj, err, N, const, scl] = ...
    pronyaaa_auto(F, dom, tol, mmax, cleanup_flag)
% pronyAAA with automated sampling of domain. 
% This is based on the aaa_auto subfunction of Chebfun's aaa function. 
% See []. 

% Flag if function has been resolved:
isResolved = 0;
% Main loop:
for n = 5:14
    % Sample points:
    Z = linspace(dom(1), dom(2), 1 + 2^n).'; Z = Z(1:end-1);
    FF = F(Z); 
    const = mean(FF); 
    FF = FF - const; 
    scl = max(abs(FF)); 
    FF = FF/scl; 
    [r, pol, res, zj, fj, wj, errvec] = pronyaaa(FF, Z, dom, tol, mmax,...
        cleanup_flag);
    
    % Test if rational interpolant is accurate:
    reltol = tol * norm(F(Z), inf);
    % on the sample: 
    err(1,1) = norm(F(Z) - (scl*r(Z)+const), inf);
    
    % on a refined grid:
    Zrefined = linspace(dom(1), dom(2), ...
        round(1.5 * (1 + 2^(n+1)))).';
    Zrefined = Zrefined(1:end-1); 
    err(2,1) = norm(F(Zrefined) - (scl*r(Zrefined)+const), inf);
    
    if ( all(err < reltol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % test on pseudo random sample points:
        xeval = rand(10,1)*(dom(2)-dom(1)) + dom(1);
        
        if ( norm(F(xeval) - (scl*r(xeval)+const), inf) < reltol )
            isResolved = 1;
            break
        end
    end
end
err = max(err); 
if ( ( isResolved == 0 ) && mmax )
    warning('RFUN:pronyaaa_auto:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end
N = length(Z); 
end 

