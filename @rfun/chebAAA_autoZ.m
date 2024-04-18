function [r, pol, res, zer, zj, fj, wj, errvec,N, const, scl] = ...
    ChebAAA_autoZ(F, dom, tol, mmax, cleanup_flag)
% Automated choice of sample set

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    FF = F(Z); 
    const = mean(FF); 
    FF = FF - const; 
    scl = max(abs(FF)); 
    FF = FF/scl;
    [r, pol, res, zer, zj, fj, wj, errvec] = chebAAA(FF, Z, 'tol', tol, ...
          'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol);

    % Test if rational approximant is accurate:
    abstol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        round(1.5 * (1 + 2^(n+1)))).';
    err(2,1) = norm(F(Zrefined) - r(Zrefined), inf);
    if ( all(err < abstol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        
        if ( norm(F(xeval) - r(xeval), inf) < abstol )
            isResolved = 1;
            break
        end
    end
end

if ( ( isResolved == 0 ) && ~mmax_flag )
    warning('AAA:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end
N = length(Z); 
end % End of AAA_AUTOZ.