function [r, pol, res, zj, fj, wj, errvec, NN] = pronyaaa(F, Z, dom, tol, mmax, cleanup_flag)
%%
% constructs a type (k-1, k) barycentric trigonometric rational interpolant to F.
% Interpolating points are chosen from Z. See [3] for more details; 
% see also rfun/constructor. 
%
% This code is based on the AAA algorithm from [1], and some pieces of the 
% code are strongly influenced by the aaa function in Chebfun [2]. 
%
%% ADD LITERATURE


cleanup_tol = 1e-10*max(F); 
% note: it isn't clear that this is the best way to handle cleanup. 
% The poles with residues of magnitude smaller than cleanup_tol are eliminated
% in the cleanup stage of the algorithm. However, there are settings where 
% these poles should be kept. 

%%
% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);
M = length(Z);NN=M;
a = dom(1); 
b = dom(2); 
%map domain to [0, 1]: 
Z = (Z-a)/(b-a); 
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
mmax = min(mmax, length(F)); 
% main AAA iteration:
for m = 1:mmax 
    if m ==1
        [~, jj] = max(abs(F)); 
    else
    % Select next support point where error is largest:
    [~, jj] = max(abs(F - R));
    end                                 % Select next support point.
    zj = [zj; Z(jj)];                   % Update support points.
    fj = [fj; F(jj)];                   % Update data values.
    J(J == jj) = [];                    % Update index vector.
    %update "Lowner" matrices:
    Ceven = [Ceven cot( 2*pi*(Z - Z(jj))./2 )];  % Next column of matrix.
    Codd =  [Codd csc( 2*pi*(Z - Z(jj))./2 )];
   
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
        % Compute weights:
        %first, we need to ensure that deg(num) = deg(den)-1; 
        % to do this, we choose must choose wj in the null 
        % space of v.' 
        if mod(m,2) %odd case
            C = Codd;
            v = fj.*exp(-1i*pi*wt_v); 
            [Q, ~] = qr(v); 
            Q = Q(:,2:end); %null space of v.' 
        else %even case
            C = Ceven;
            v = fj; 
            [Q, ~] = qr(v); 
            Q = Q(:,2:end); 
        end
      
        %compute weights:
        Sf = diag(fj);                   
        A = (SF*C - C*Sf);                
        [~, ~, V] = svd(A(J,:)*Q, 0);         
        wj = V(:,end);     
        wj = Q*wj;          
     
        % Evaluate rational on Z:
        N = C*(wj.*fj);                    
        D = C*wj;                           
        R = F;
        R(J) = N(J)./D(J);
    end
        % Error in the sample points:
        err = norm(F - R, inf);
        errvec = [errvec; err];
    
    % Check if converged (also, we want even # of sampled points):
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
 

% Construct function handle to evaluate r:
r = @(zz) reval(zz, zj, fj, wj);

% Compute poles, residues:
[pol, res] = pr(zj, fj, wj, dom);

if (cleanup_flag)
    % Remove spurious poles:
    [r, pol, res, zj, fj, wj] = cleanup(r, pol, res, zj, fj, wj,...
        Z, F, dom,cleanup_tol);
end

% sort poles into conjugate pairs (if possible):

pol1 = pol(imag(pol)>0); 
pol2 = pol(imag(pol)<0); 
if length(pol1)==length(pol2) && max(abs(imag(pol))) < 1e-12
 pol1 = sort(pol1,'ComparisonMethod',  'real'); 
 pol2 = sort(pol2, 'ComparisonMethod', 'real'); 
 pol(1:2:end) = pol1; 
 pol(2:2:end) = pol2; 
end

%adjust nodes: 
zj = zj*(b-a) + a; 

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate rational function in barycentric form.

function r = reval(zz, zj, fj, wj)
% Evaluate rational function in barycentric form.
 m = length(zj); 
 zv = zz(:);  

if mod(m,2)
    CC = bsxfun(@(vals, params) csc(2*pi*(vals - params)./2), zv, zj.');
else
    CC = bsxfun(@(vals, params) cot(2*pi*(vals - params)./2), zv, zj.'); 
end

r = (CC*(wj.*fj))./(CC*wj);             % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL().


%% Compute poles, residues of the rational.

function [pol, res] = pr(zj, fj, wj, dom)
m = length(wj);
a = dom(1); 
b = dom(2); 
% Compute poles via generalized eigenvalue problem:
B = diag([ones(m, 1); 0]); 
B(1:m, end) = -1i*wj; 
E = diag([exp(1i*2*pi*zj);0]); 
E(1:m,end) = 1i*wj.*exp(1i*2*pi*zj); 
E(end, 1:m) = ones(1, m); 
pol = eig(E, B); 
pol = log(pol)/(2i*pi);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

%correct poles so that they are on the interval [a, b]: 
pol = a + pol*(b-a); 
pol = a + mod(real(pol), b-a) + 1i*imag(pol);
zj = a + zj*(b-a); 
%%
% Compute residues via l'Hopital's rule: res f(c) = n(c)/d'(c)
CC1 = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), pol, zj.'); 
CC2 = bsxfun(@(vals, params) (csc(2*pi/(b-a)*(vals - params)./2)).^2, pol, zj.'); 

nc = CC1*(wj.*fj);
nd = -pi/(b-a)*CC2*wj; 

res = nc./nd; 

end 


%% Cleanup:Remove spurious poles (including any real-valued poles.)

function [r, pol, res, z, f, w] = cleanup(r, pol, res, z, f, w, Z, F, dom, tol)
% this routine tries to eliminate (i) poles with small residues and (ii)
% any real-valued poles.
a = dom(1); 
b = dom(2); 
ii = find(abs(res) < 2*tol* norm(F, inf) | abs(imag(pol)) < 1e-11);
ni = length(ii);
if ~isempty(ii) && mod(ni, 2)
    %find the next smallest residue: 
    [~, idx] = sort(abs(res)); 
    ii = [ii; idx(length(ii)+1)];  
    ni = length(ii); 
end
if ( ni == 0 ) || length(res)==2
    return
else
    fprintf('%d Froissart doublets.\n', ni)
end

% For each spurious pole, we  find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
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

M = length(Z);
% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
C = bsxfun(@(vals, params) cot(2*pi*(vals - params)./2), Z, z.'); 
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[Q, ~] = qr(f); 
Q = Q(:,2:end); 
[~, ~, V] = svd(A*Q, 0);
w = V(:,end);
 w = Q*w; 
% Build function handle and compute poles, residues:
r = @(zz) reval(zz, z, f, w);
[pol, res] = pr( z, f, w, dom);
end % End of CLEANUP().






