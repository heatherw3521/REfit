function [r, pol, res, zer, zj, fj, wj, errvec, N] = pronyaaa_poles(F,P, Z, dom, tol, mmax, cleanup_flag)
%   Computes a type (k, k-1) rational trigonometric rational approximation
%   that has P as a subset of suggested poles.
%   
%   See [], and also efun/ift. 

%% NOTE: IS THIS FUNCTION OUTDATED/NOT USEFUL?

%% ADD LITERATURE


cleanup_tol = 1e-10; 
lp = length(P);

% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);
M = length(Z);

% Relative tolerance:
reltol = tol * norm(F, inf);
% Left scaling matrix:
SF = spdiags(F, 0, M, M);

% Initialization for AAA iteration:
J = 1:M;
zj = [];
fj = [];
Ceven = [];
Codd = []; 
errvec = [];
R = mean(F);
wt_v = []; %odd weight null space adjustment

% AAA iteration:
for m = 1:mmax    
    % Select next support point where error is largest:
    [~, jj] = max(abs(F - R));          % Select next support point.
    zj = [zj; Z(jj)];                   % Update support points.
    fj = [fj; F(jj)];                   % Update data values.
    J(J == jj) = [];                    % Update index vector.
    %update "Lowner" matrices:
    Ceven = [Ceven cot( 2*pi/(dom(2)-dom(1))*(Z - Z(jj))./2 )];  % Next column of "Cauchy" matrices.
    Codd =  [Codd csc( 2*pi/(dom(2)-dom(1))*(Z - Z(jj))./2 )];
   
    %update vector for controlling degree of num: 
    if m == 1 
        wt_v = []; 
    elseif m==2
        wt_v = [zj(2); zj(1)]; 
    else
        wt_v = wt_v + zj(end); 
        wt_v = [wt_v; sum(zj(1:end-1))]; %add next col
    end
    
    if m == 1 %the first (m-1, m) rational function is zero:
            % just compute error manually. 
        R = F; 
        R(J) = 0; 
    else
        % Compute weights (we use poles for this):
        if mod(m,2) %odd case
            C = Codd;
            v = fj.*exp(-1i*pi/(dom(2)-dom(1))*wt_v); 
            VV = bsxfun(@(vals, params) csc(2*pi/(dom(2)-dom(1))*(vals - params)./2), P, zj.'); 
            VV = [VV; v.']; 
        else %even case
            C = Ceven;
            v = fj; 
            VV = bsxfun(@(vals, params) cot(2*pi/(dom(2)-dom(1))*(vals - params)./2), P, zj.'); 
            VV = [VV; v.']; 
        end
        if m <= lp %if we have less poles than points, we solve w/VV directly:
            [~, ~, V] = svd(VV, 0); 
            wj = V(:,end); 
        else %solve the least squares problem but restrict solutions to null space of VV
        [Q, ~]= qr(VV'); 
        Q = Q(:,lp+1:end); 
        Sf = diag(fj);                      % Right scaling matrix.
        A = (SF*C - C*Sf);                % Loewner matrix
        [~, ~, V] = svd(A(J,:)*Q, 0);         % Reduced SVD.
        wj = V(:,end);      % weight vector in null space of v = min sing vector
        wj = Q*wj;          % compute weight vector
        end
        % Rational approximant on Z:
        N = C*(wj.*fj);                     % Numerator
        D = C*wj;                           % Denominator
        R = F;
        R(J) = N(J)./D(J);
    end
        % Error in the sample points:
        err = norm(F - R, inf);
        %errvec = [errvec; err];   
    % Check if converged: (also, we want even # of sampled points)
    if ( err <= reltol && mod(m,2)==0)
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
N = length(F); %length of domain sample. 
% Construct function handle:
r = @(zz) reval(zz, zj, fj, wj, dom);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(zj, fj, wj, dom);

if (cleanup_flag)
    % Remove Froissart doublets:
    [r, pol, res, zer, zj, fj, wj] = cleanup(r, pol, res, zer, zj, fj, wj,...
        Z, F, dom,cleanup_tol);
end

% sort poles into conjugate pairs (if possible):

pol1 = pol(imag(pol)>0); 
pol2 = pol(imag(pol)<0); 
%pol3 = pol(imag(pol)==0); 
% if ~(isempty(pol3))  %this pole shouldn't exist!
%     error('ran into real-valued pole')
%     %to do: add a 'remove pole' command
% else 
if length(pol1)==length(pol2) && max(abs(imag(pol))) < 1e-12
 pol1 = sort(pol1,'ComparisonMethod',  'real'); 
 pol2 = sort(pol2, 'ComparisonMethod', 'real'); 
 pol(1:2:end) = pol1; 
 pol(2:2:end) = pol2; 
end
 %to do: throw poles that are not conj pairs on the end

end 



%% parse Inputs:

function [F, Z, M, dom, tol, mmax, cleanup_flag, needZ, mmax_flag] = ...
    parseInputs(F, varargin)
% Input parsing for AAA.

% Check if F is empty:
if ( isempty(F) )
    error('RFUN:aaa:emptyF', 'No function given.')
elseif ( isa(F, 'chebfun') )
    if ( size(F, 2) ~= 1 )
        error('RFUN:aaa:nColF', 'Input chebfun must have one column.')
    end
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('RFUN:aaa:emptyZ', ...
            'If sample set is provided, it must be nonempty.')
    end
    varargin(1) = [];
end

% Set defaults for other parameters:
tol = 1e-10;        % Relative tolerance.
mmax = 200;         % Maximum number of terms.
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [0, 1];
end
cleanup_flag = 1;   % Cleanup on.
mmax_flag = 0;

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            mmax = varargin{2};
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('RFUN:aaa:dom', ...
                    ['Given domain does not match the domain of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
        if ( strncmpi(varargin{2}, 'off', 3) || ( varargin{2} == 0 ) )
            cleanup_flag = 0;
        end
        varargin([1, 2]) = [];
        
    else
        error('RFUN:aaa_trig:UnknownArg', 'Argument unknown.')
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
            error('RFUN:aaa:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    elseif ( ischar(F) )
        % F is given as a string input. Convert it to a function handle.
        F = str2op(vectorize(F));
        F = F(Z);
    else
        error('RFUN:aaa:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    needZ = 1;
    Z = [];
    M = length(Z);
end

end % End of PARSEINPUT().


%% Evaluate rational function in barycentric form.

function r = reval(zz, zj, fj, wj, dom)
% Evaluate rational function in barycentric form.
 m = length(zj); 
 b = dom(2); 
 a = dom(1); 
 zv = zz(:);  

if mod(m,2)
    CC = bsxfun(@(vals, params) csc(2*pi/(b-a)*(vals - params)./2), zv, zj.');
else
    CC = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), zv, zj.'); 
end

r = real((CC*(wj.*fj))./(CC*wj));             % vector of values

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

end % End of REVAL().


%% Compute poles, residues and zeros.

function [pol, res, zer] = prz(r, zj, fj, wj, dom)
% Compute poles, residues, and zeros of rational function in barycentric form.
m = length(wj);
a = dom(1); 
b = dom(2); 
% Compute poles via generalized eigenvalue problem:
B = diag([ones(m, 1); 0]); 
B(1:m, end) = -1i*wj; 
E = diag([exp(1i*2*pi/(dom(2) - dom(1))*zj);0]); 
E(1:m,end) = 1i*wj.*exp(1i*2*pi/(dom(2) - dom(1))*zj); 
E(end, 1:m) = ones(1, m); 
pol = eig(E, B); 
pol = log(pol)/(2i*pi/(dom(2) - dom(1)));
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

%%


% sort poles into conjugate pairs:
%deal with poles close to axis (probably spurious): 

% pol1 = pol(imag(pol)>1e-14); 
% pol2 = pol(imag(pol)<1e-14); 
% pol3 = pol(imag(pol)==0); 
% if ~(isempty(pol3))  %this pole shouldn't exist!
%     error('ran into real-valued pole')
%     %to do: add a 'remove pole' command
% else    
% pol1 = sort(pol1,'ComparisonMethod',  'real'); 
% pol2 = sort(pol2, 'ComparisonMethod', 'real'); 
% pol(1:2:end) = pol1; 
% pol(2:2:end) = pol2; 
% end

%should we average and make these truly complex conjugates?
%%
% Compute residues via discretized Cauchy integral:
%dz = 1e-5*exp(2i*pi*(1:4)/4);
%res = r(bsxfun(@plus, pol, dz))*dz.'/4;


% Compute residues via l'Hopital's rule: res f(c) = n(c)/d'(c)
CC1 = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), pol, zj.'); 
CC2 = bsxfun(@(vals, params) (csc(2*pi/(b-a)*(vals - params)./2)).^2, pol, zj.'); 

nc = CC1*(wj.*fj);
nd = -pi/(b-a)*CC2*wj; 

res = nc./nd; 

% Compute zeros via generalized eigenvalue problem:
E(end, 1:m) = fj;
zer = eig(E, B);
zer = log(zer)/(2*pi/(dom(2) - dom(1))*1i);
% Remove zeros of numerator at infinity:
%remove zeros at z = \pm infinity and z = 0:
zer = zer(abs(zer)> 1e-10); 
zer = zer((abs(zer)<=1e10)); 
zer = zer(~isinf(zer));
if length(zer) > m-2 %check that the zero eig is dropped
    zer = zer(abs(zer) > 1e-6); 
end


end % End of PRZ().


%% Cleanup

function [r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F, dom, tol)
% Remove spurious pole-zero pairs.

% Note: we want support points removed in pairs, since we should 
% always have an even # of points to preserve parity of rational function
 
a = dom(1); 
b = dom(2); 
%res_both = abs(res(1:2:end)+res(2:2:end)); %just look at one in the pair; 
% Find negligible residues:
%ii = find(res_both < 2*tol* norm(F, inf));
ii = find(abs(res) < 2*tol* norm(F, inf));
%idx = 1:2:length(res); 
%ii = [idx(ii); idx(ii)+1]; %pick out the pairs; 
%ii = sort(ii); 
ni = length(ii);
if ~isempty(ii) && mod(ni, 2)
    warning('RFUN:constructor:aaa:Odd number of Froisssart doublets detected.')
    %find the next smallest residue: 
    [~, idx] = sort(abs(res)); 
    ii = [ii; idx(length(ii)+1)];  
    ni = length(ii); 
end
if ( ni == 0 ) || length(res)==2
    % Nothing to do.
    return
else
    fprintf('%d Froissart doublets.\n', ni)
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
% Remove support points z from sample set:
    %F(Z == z(jj)) = [];
    %Z(Z == z(jj)) = [];
    % Remove support point(s):
    z(jj) = [];
    f(jj) = [];
    w(jj) = []; %match removals
end

%now remove the rest of the support points from sample: 
for j = 1:length(z)
    F(Z == z(j)) = [];
    Z(Z == z(j)) = [];
end
    

%if ni is odd, we need to remove one more support point to make
% an even set: choose the smallest residue:
% (note, we need to make this more sophisticated ?)

%if mod(ni, 2)
   % mj = find(rr == min(abs(res)),1); 
   % z(mj) = []; 
   % f(mj) = []; 
   % F(Z == z(mj)) = []; 
   % Z(Z == z(mj)) = []; 
%end

m = length(z);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
%C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
C = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), Z, z.'); 
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[Q, ~] = qr(f); 
Q = Q(:,2:end); 
[~, ~, V] = svd(A*Q, 0);
w = V(:,end);
 w = Q*w; 
% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w, dom);
[pol, res, zer] = prz(r, z, f, w, dom);
%call again to check for more doublets:
%[r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F, dom, tol);
end % End of CLEANUP().


%% Automated choice of sample set

function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, mmax_flag)
%

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    fprintf('Sampling \n')
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa_trig(F, Z, 'tol', tol, ...
        'mmax', mmax, 'cleanup', cleanup_flag);
    
    % Test if rational approximant is accurate:
    reltol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        round(1.5 * (1 + 2^(n+1)))).';
    err(2,1) = norm(F(Zrefined) - r(Zrefined), inf);
    
    if ( all(err < reltol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        
        if ( norm(F(xeval) - r(xeval), inf) < reltol )
            isResolved = 1;
            break
        end
    end
end

if ( ( isResolved == 0 ) && ~mmax_flag )
    warning('RFUN:aaa:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end

end % End of AAA_AUTOZ().

function op = str2op(op)
    % Convert string inputs to either numeric format or function_handles.
    sop = str2num(op);
    if ( ~isempty(sop) )
        op = sop;
    else
        depVar = symvar(op);
        if ( numel(depVar) ~= 1 )
            error('CHEBFUN:CHEBFUN:str2op:indepvars', ...
             'Incorrect number of independent variables in string input.');
        end
        op = eval(['@(' depVar{:} ')', op]);
    end
end % End of STR2OP().



