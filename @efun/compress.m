function s = compress(S, varargin)
% When S is an efun with a suboptimal number of terms, 
% this function returns s, a shorter exponential sum. 
%
% tolerance can be specified with compress(S, tol). 
% Default tol = S.tol. 
%

% See also: efun/plus, efun/times, rfun/convolve.

%%
%%empty case:
if isempty(S) 
    %pass to constructor 
    s = efun(); 
    return
end

m = length(S);
%length 1 case: 
if (m==1) 
    s = S; 
    return
end

tol = S.tol; 
if ~isempty(varargin)   
    tol = varargin{1};   
end

%% rectangular Prony-like method. 
happy = false; 
cp = 0;
const = S.const; 
S.const = 0; 
scl = S.scl; 
S.scl = 1; 
D = S.res; % we want accuracy up to bandlimit D
s = S;  
M = 4*m+1; %initial sample size. 
while ~happy
    cp = [cp; feval(S,(length(cp):M).')];
    %h = efun(feval(S,(0:M).', 'coeffs', 'length', M, 'tol', tol, 'rectangle');
    [r,ss,w] = prony_compress(cp, m, tol); 
    s.exp = r; 
    s.weights = w; 
    s.const = 0; 
   % test accuracy; 
    rnd = randi([0 D], 1,50).'; 
    err = abs(feval(s,[(0:length(cp)-1).';rnd])-feval(S,[(0:length(cp)-1).';rnd]));
    if max(err) < tol 
        happy = true;
        break;
    end
    if M > 2*D 
        break; 
    end 
    M = 2*M; %double # samples; try again
end
    %% 
    s.sv = ss; 
    s.tol = tol; 
    s.const = const; 
    s.scl = scl; 
end
    






