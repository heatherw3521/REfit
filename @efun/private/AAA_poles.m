function [r, pol, res, zer, zj, fj, wj] = AAA_poles(F,Z,P,xp,fj, varargin)
%   A version of AAA where the poles are specified. 
%   See Wilber, Damle, Townsend, 2021 for details. 
%%

% TO DO: add a test to find real-valued poles and then resample if needed. 
%compute the barycentric weights
dom = [0, 1]; 
 
lp = length(P);
m= length(xp);

% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);

%remove points xp from sample
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);

if any(ismember(Z, xp))
    [~, jj]= ismember(xp, Z); 
    F(jj(jj>0)) = []; 
    Z(jj(jj>0)) = []; 
end
zj = xp; 
M = length(F); 

% Left scaling matrix:
SF = spdiags(F, 0, M, M);
C = bsxfun(@(vals, params) cot( 2*pi/(dom(2)-dom(1))*(vals - params)./2 ), Z, zj.'); 
v = fj; 
VV = bsxfun(@(vals, params) cot(2*pi/(dom(2)-dom(1))*(vals - params)./2), P, zj.'); 
VV = [VV; v.'];
%%

[Q, ~]= qr(VV'); 
Q = Q(:,lp:end);  
Sf = diag(fj);                      % Right scaling matrix.
A = (SF*C - C*Sf);                % Loewner matrix
[~, ~, V] = svd(A*Q, 0);         % Reduced SVD.
wj = V(:,end);      % weight vector in null space of v = min sing vector
wj = Q*wj;
 
% Remove support points with zero weights: 
I = find(wj == 0);
if ~isempty(I) && mod(length(I), 2)
    [~, fm] = min(abs(wj)); %for odd # of zeros, we need to elimimate one more 
    I = [I(:); fm]; 
end

% Compute poles, residues and zeros:
[pol, res, zer] = prz( zj, fj, wj, dom);

%remove smallest residues: 
 
a = dom(1); 
b = dom(2); 

[~,ii] = sort(abs(res));

ni = m-lp;
if ( ni==0) || length(res)==2 
    % Nothing to do.
    r = @(zz) reval(zz, zj, fj, wj, dom);
    [pol, res, zer] = prz( zj, fj, wj, dom);
    return
%else
    %fprintf('%d Froissart doublets.\n', ni)
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(zj-pol(ii(j)));
    jj = find(azp == min(azp),1);
% Remove support points z from sample set:
    % Remove support point(s):
    zj(jj) = [];
    fj(jj) = [];
    wj(jj) = []; %match removals
end

%now remove the rest of the support points from sample: 
for j = 1:length(zj)
    F(Z == zj(j)) = [];
    Z(Z == zj(j)) = [];
end

m = length(zj);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(fj);
C = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), Z, zj.'); 
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[Q, ~] = qr(fj); 
Q = Q(:,2:end); 
[~, ~, V] = svd(A*Q, 0);
wj = V(:,end);
wj = Q*wj; 
% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, zj, fj, wj, dom);
[pol, res, zer] = prz( zj, fj, wj, dom);
end

%% Compute poles, residues and zeros.

function [pol, res, zer] = prz( zj, fj, wj, dom)
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
%remove zeros at z = \pm infinity and z = 0:
zer = zer(abs(zer)> 1e-10); 
zer = zer((abs(zer)<=1e10)); 
zer = zer(~isinf(zer));
if length(zer) > m-2 %check that the zero eig is dropped
    zer = zer(abs(zer) > 1e-6); 
end

%put in terms of x
zer = log(zer)/(2*pi/(dom(2) - dom(1))*1i);
%move to interval: 

end % End of PRZ().

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
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL().