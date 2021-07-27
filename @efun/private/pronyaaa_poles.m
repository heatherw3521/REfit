function [pol, res, zj, fj, wj] = pronyaaa_poles(F,Z,P,dom, pp)
% this function performs a version of AAA designed to 
% recover a trigonometric rational R from samples F at locations Z, where
% the poles R are given by P. 
%
% See [1] for details. Also see: efun/ift, rfun/private/pronyaaa


%%
% STEP 1: choose candidate support points using CPQR pivoting
P = P(:); %poles of R
Z = Z(:); %sample locations

a = dom(1); 
b = dom(2); 
%shift samples and poles to [0, 1): 
Z = (Z-a)/(b-a); 
P = (P-a)/(b-a); 
F = F(:); %samples of R
%pp = oversampling parameter
lp = length(P);
m = lp+pp; 
 
% Remove any infinite or NaN function values (avoid SVD failures):
toKeep = ~isinf(F);
%remove points xp from sample
F = F(toKeep); Z = Z(toKeep);
toKeep = ~isnan(F);
F = F(toKeep); Z = Z(toKeep);
  
% set up matrix (denominator of Barycentric form): 
V = cot(pi*(P - Z.'));
V = [V; F.']; 
[~, ~, idx] = qr(V, 'vector');
%
% if the gap is large, we can solve directly for the weights using V: 
% rdiag = diag(R)/R(1,1); 
% if rdiag(end)/rdiag(end) > 1e-2 && rdiag(end)<1e-11
%     [~, ~, W] = svd(V(:,idx(1:lp))); 
%     wj = W(:, end);  
%     zj = Z(idx(1:lp)); 
%     fj = F(idx(1:lp)); 
%     [pol, res] = pr(zj, fj, wj, dom);
%     zj = zj*(b-a)+a; %correct nodes to match domain
%     return
% end
zj = Z(idx(1:m)); %candidate support points 
Zo = Z; 
Z(idx(1:m)) = []; 
fj = F(idx(1:m)); %values at candidate points 
Fo = F; 
F(idx(1:m)) = []; 
M = length(F); 
  
%%
% STEP 2: solve a modified linear least squares problem. 
  
%Left scaling matrix
SF = spdiags(F, 0, M, M);
C = bsxfun(@(vals, params) cot( pi*(vals - params) ), Z, zj.'); 
VV = bsxfun(@(vals, params) cot(pi*(vals - params)), P, zj.'); 
VV = [VV; fj.'];

%projection:
[Q, ~]= qr(VV'); 
Q = Q(:,lp:end); 

Sf = diag(fj);                 
A = (SF*C - C*Sf);               
[~, ~, V] = svd(A*Q, 0);        
wj = V(:,end);      % weight vector
wj = Q*wj;
  
% Compute poles and residues:
[pol, res] = pr(zj, fj, wj, dom);

%now find smallest residues in order to subselect 
% support points from zj:
a = dom(1); 
b = dom(2); 
[~,ii] = sort(abs(res));

ni = m-lp;
if (ni==0) || length(res)==2 
    % No correction to make.
    [pol, res] = pr(zj, fj, wj, dom);
    return
end

% now find and remove closest support point:
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
%for j = 1:length(zj)
   % F(Z == zj(j)) = [];
   % Z(Z == zj(j)) = [];
%end

%recompute weights:
%M = length(Z);
SF = spdiags(F, 0, M, M);
Sf = diag(fj);
C = bsxfun(@(vals, params) cot(2*pi*(vals - params)./2), Z, zj.'); 
A = SF*C - C*Sf;                    

[Q, ~] = qr(fj); 
Q = Q(:,2:end); 
[~, ~, V] = svd(A*Q, 0);
wj = V(:,end);
wj = Q*wj; 
% Build function handle and re-compute poles, residues:
[pol, res] = pr(zj, fj, wj, dom);

% check for real-valued poles: 
if any(abs(imag(pol)) < 1e-10)
    pp = pp + 2; % try again with a higher oversampling parameter
    [pol, res] = cleanup(Fo, Zo,P, idx, dom,pp); 
end
zj = zj*(b-a)+a; %correct nodes to match domain
end

%%%%%%%%%%%%%%%%%%%%%%END MAIN%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pol, res] = cleanup(F, Z, P,idx, dom, pp)

lp = length(P);
m = lp+pp; 

zj = Z(idx(1:m)); %candidate support points  
Z(idx(1:m)) = []; 
fj = F(idx(1:m)); %values at candidate points  
F(idx(1:m)) = []; 
M = length(F); 
  
%%
% STEP 2: solve a modified linear least squares problem. 
  
%Left scaling matrix
SF = spdiags(F, 0, M, M);
C = bsxfun(@(vals, params) cot( pi*(vals - params)), Z, zj.'); 
%projection: (TO DO: create update Q to make this faster, since we just need to append two more orthog cols)
VV = bsxfun(@(vals, params) cot(pi*(vals - params)), P, zj.'); 
VV = [VV; fj.'];
[Q, ~]= qr(VV'); 
Q = Q(:,lp:end); 

Sf = diag(fj);                  % Right scaling matrix.
A = (SF*C - C*Sf);               
[~, ~, V] = svd(A*Q, 0);        
wj = V(:,end);      % weight vector
wj = Q*wj;
  
% Compute poles and residues:
[pol, res] = pr(zj, fj, wj, dom);

%now remove smallest residues in order to subselect 
% support points from zj:
[~,ii] = sort(abs(res));

ni = m-lp;
if (ni==0) || length(res)==2 
    % No correction to make.
    [pol, res] = pr(zj, fj, wj, dom);
    return
end
% now find and remove closest support point:

for j = 1:ni
    azp = abs(zj-pol(ii(j)));
    jj = find(azp == min(azp),1);
% Remove support points z from sample set:
    % Remove support point(s):
    zj(jj) = [];
    fj(jj) = [];
    wj(jj) = []; %match removals
end

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(fj);
C = bsxfun(@(vals, params) cot(pi*(vals - params)), Z, zj.'); 
A = SF*C - C*Sf;                    

% Solve least-squares problem to obtain weights:
[Q, ~] = qr(fj); 
Q = Q(:,2:end); 
[~, ~, V] = svd(A*Q, 0);
wj = V(:,end);
wj = Q*wj; 
% Build function handle and compute poles, residues and zeros:
[pol, res] = pr(zj, fj, wj, dom);

end

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

%correct poles and nodes so that they are on the interval [a, b]: 
pol = mod(real(pol), 1) + 1i*imag(pol); 
pol = pol*(b-a)+a;
zj = zj*(b-a)+a;
%%
% Compute residues via l'Hopital's rule: res f(c) = n(c)/d'(c)
CC1 = bsxfun(@(vals, params) cot(2*pi/(b-a)*(vals - params)./2), pol, zj.'); 
CC2 = bsxfun(@(vals, params) (csc(2*pi/(b-a)*(vals - params)./2)).^2, pol, zj.'); 

nc = CC1*(wj.*fj);
nd = -pi/(b-a)*CC2*wj; 

res = nc./nd; 

end 

  
  

  
      

