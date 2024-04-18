function [r, pol, res, zer, zj, fj, wj, errvec,N] = chebAAA(F, varargin)
% aaa rational approximation of data F on set Z
% [r,pol,res,zer,z,f,w,errvec] =chebAAA(F,Z,tol,mmax)
%
% Input: F = vector of data values, or a function handle
% Z = vector of sample points
% tol = relative tolerance tol, set to 1e-13 if omitted
% mmax: max type is (mmax-1,mmax), set to 100 if omitted
%
% Output: r = AAA approximant to F (function handle)
% pol,res,zer = vectors of poles, residues, zeros
% z,f,w = vectors of support pts, function values, weights
% errvec = vector of errors at each step

%% initialization
%phrase input
[F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, needZ, mmax_flag, ...
    nlawson, degree_flag, degree] = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, ...
            nlawson, degree_flag, degree);
    return
end
% Remove infinite or NaN function values and repeated entries:
toKeep = ~isinf(F) & ~isnan(F);
F = F(toKeep); Z = Z(toKeep);
[Z, uni] = unique(Z,'stable'); F = F(uni);

M = length(Z);N=M;
abstol = tol*norm(F, inf);                 % Absolute tolerance
J = (1:M)';
zj = []; fj = []; C = []; A = [];
errvec = [];
R = mean(F)*ones(size(J));
SF = spdiags(F,0,M,M); % left scaling matrix

J = (1:M)';
z = [];
f = [];
C = []; 
errvec = []; 
R = mean(F)*ones(size(J));
%% main AAA loop
for m = 1:mmax
    [~,jj] = max(abs(F(J) - R(J)));       % Select next support point
    zj = [zj; Z(J(jj))];                   % Update support points
    fj = [fj; F(J(jj))];                   % Update data values
    C = [C 1./(Z - Z(J(jj)))];             % Next column of Cauchy matrix
    J(jj) = [];                            % Update index vector
    A = [A, (F-fj(end)).*C(:,end)];        % Update Loewner matrix
    
    if m == 1 %the first (m-1, m) rational function is zero:
            % just compute error manually. 
        R = F; 
        R(J) = 0; 
    else
        % Qr iteration to make sure the function is (n-1,n)
        [Q, ~] = qr(fj);  %QR step
        Q = Q(:,2:end);   
        [~, ~, V] = svd(A(J,:)*Q, 0);   %least square on |CQz|      
        wj = V(:,end);       % find best z
        wj = Q*wj;           %compute weight w = Qz
        
        N = C*(wj.*fj);                    
        D = C*wj;                           
        R = F;
        R(J) = N(J)./D(J);
    end
        % Error in the sample points:
        err = norm(F - R, inf);
        errvec = [errvec; err];
    
         % Check if converged (also, we want even # of sampled points):
    maxerr = norm(F - R, inf);
    errvec = [errvec; maxerr];
    if ( maxerr <= abstol )
        break
    end
end
    

% Remove support points with zero weight:
    I = find(wj == 0);
    if ~isempty(I) && mod(length(I), 2)
        [~, fm] = min(abs(wj)); %for odd # of zeros, we need to elimimate one more 
        I = [I(:); fm]; 
    end
    
    zj(I) = [];
    wj(I) = [];
    fj(I) = [];

    r = @(zz) reval(zz, zj, fj, wj);
    [pol, res, zer] = prz(zj, fj, wj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   AAA_AUTOZ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, cleanup_tol, mmax_flag, ...
               nlawson, degree_flag, degree)
% Automated choice of sample set

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    if degree_flag
       [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
          'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, ...
          'lawson', nlawson, 'degree', degree);
    else
       [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
          'mmax', mmax, 'cleanup', cleanup_flag, 'cleanuptol', cleanup_tol, ...
          'lawson', nlawson);
    end
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

end % End of AAA_AUTOZ.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PARSEINPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, Z, M, dom, tol, mmax, cleanup_flag, cleanup_tol, ...
    needZ, mmax_flag, nlawson, degree_flag, degree] = parseInputs(F, varargin)

% Check if F is empty:
if ( isempty(F) )
    error('AAA:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('AAA:nColF', 'Input chebfun must have one column.')
    end
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('AAA:emptyZ', ...
            'If sample set is provided, it must be nonempty.')
    end
    varargin(1) = [];
end

% Set defaults for other parameters:
tol = 1e-13;                   % Relative tolerance
mmax = 100;                    % Maximum number of terms
degree = NaN;                  % Specified degree
cleanup_tol = 1e-13;           % Cleanup tolerance
nlawson = Inf;                 % Number of Lawson steps (Inf means adaptive)
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end
cleanup_flag = 1;              % Cleanup on
mmax_flag = 0;                 % Checks if mmax manually specified
degree_flag = 0;               % Checks if degree specified
cleanup_set = 0;               % Checks if cleanup_tol manually specified
while ( ~isempty(varargin) )   % Check if parameters have been provided
    if ( strncmpi(varargin{1}, 'tol', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol = varargin{2};
            if ( ~cleanup_set && tol > 0 ) % If not set, set cleanup_tol to tol
              cleanup_tol = tol;
            end
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'degree', 6) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2}+1 )
                error('AAA:degmmaxmismatch', ' mmax must equal degree+1.')
            end            
            degree = varargin{2};
            mmax = degree + 1;
            mmax_flag = 1; 
            degree_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )            
            if ( mmax_flag == 1 ) && ( mmax ~= varargin{2})                
                error('AAA:degmmaxmismatch', ' mmax must equal degree+1.')
            end
            mmax = varargin{2};
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'lawson', 6) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            nlawson = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('AAA:dom', ...
                    ['Given domain does not match that of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif ( strncmpi(varargin{1}, 'cleanuptol', 10) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
          cleanup_tol = varargin{2};
          cleanup_set = 1;
        end
        varargin([1, 2]) = [];

    elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
        if ( strncmpi(varargin{2}, 'off', 3) || ( varargin{2} == 0 ) )
            cleanup_flag = 0;
        elseif ( varargin{2} == 2 )     % Alternative cleanup
            cleanup_flag = 2;
        end
        varargin([1, 2]) = [];
        
    else
        error('AAA:UnknownArg', 'Argument unknown.')
    end
end

% Deal with Z and F:
if ( ~exist('Z', 'var') && isfloat(F) )
    % F is given as data values, pick same number of sample points:
    Z = linspace(dom(1), dom(2), length(F)).';
end

if ( exist('Z', 'var') )
    % Z is given:
    needZ = 0;
    
    % Work with column vector:
    Z = Z(:);
    M = length(Z);
    
    % Function values:
    if ( isa(F, 'function_handle') || isa(F, 'chebfun') )
        % Sample F on Z:
        F = F(Z);
    elseif ( isnumeric(F) )
        % Work with column vector and check that it has correct length.
        F = F(:);
        if ( length(F) ~= M )
            error('AAA:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    elseif ( ischar(F) )
        % F is given as a string input. Convert it to a function handle.
        F = inline(vectorize(F));
        F = F(Z);
    else
        error('AAA:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    needZ = 1;
    Z = [];
    M = length(Z);
end

if ( ~degree_flag && (nlawson == Inf) )
    nlawson = 0;               
end

end % End of PARSEINPUTS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PRZ   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pol, res, zer] = prz(zj, fj, wj)
%   Compute poles, residues, and zeros of rational fun in barycentric form.

% Compute poles via generalized eigenvalue problem:
m = length(wj);
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
pol = pol(~isinf(pol));

% Compute residues via formula for res of quotient of analytic functions:
N = @(t) (1./(t-zj.')) * (fj.*wj);
Ddiff = @(t) -((1./(t-zj.')).^2) * wj;
res = N(pol)./Ddiff(pol);

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
zer = zer(~isinf(zer));

end % End of PRZ.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REVAL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = reval(zz, zj, fj, wj)
%   Construct function handle to evaluate rational function in barycentric form.

zv = zz(:);                         % vectorize zz if necessary
CC = 1./(zv-zj.');                  % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);         % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL.