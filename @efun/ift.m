function r = ift(s,varargin)
% The inverse Fourier transform of s.  
%
% r = ift(s) returns an rfun representing a barycentric trigonometric rational 
% function corresponding to the inverse Fourier transform of s (See [2]). 
%
% ift(s, 'pr') returns a handle to the pole-residue form of r(x).
% ift(s, 'zt') returns a handle for the pole-residue form of 
% r(z), where z = exp(2*pi*1i*x/diff(dom(s.domain)). 
%
% For rfuns:
% ift(s, 'cpqr')  uses CPQR pivoting to do AAA (default) (See [2]).
%
% ift(s, 'aaa') applies standard AAA (i.e., greedy pivoting) (See [1]).
%
% ift(s, locs, 'type') applies AAA on the grid given by locs, 
% with pivoting style given by 'type'. 
%
% See also rfun/ft. 
%%

%% ADD LITERATURE

if isempty(s)
    r = []; 
    return
end

%% 
% set defaults:
type = 'cpqr';  
tol = s.tol; 
fft_flag = 0; 
pts = []; 
nums = []; 

%parse input: 
for j = 1:nargin-1
    lk = varargin{j};
    if isa(lk, 'double')
        if all(size(lk)==1) %input is a sample size
            nums = lk + mod(lk, 2); 
        else
            pts = lk(:); %input is a grid to sample on. 
        end
    elseif strcmpi(lk, 'fft') || strcmpi(lk, 'ifft') %use fft to sample
        fft_flag = 1; 
    elseif strcmpi(lk, 'aaa') || strcmpi(lk, 'cpqr') || strcmpi(lk, 'polres') || strcmpi(lk, 'zt')
        %%input tells us the type of ift
        type = lk; 
    else
        error('efun:ift:cannot parse input')
    end
end
%%
% deal with types asking for function handles: 
if (strcmpi(type, 'polres')  || strcmpi(type,'zt') )
    dom = s.domain;
    pol = poles(s, 'zt'); pol = pol(abs(pol)<1); 
    pol = [pol ; 1./conj(pol)];
    w = s.weights; 
    const = s.const; 
    scl = s.scl; 
    res = [w.*exp(s.exp); -conj(w).*exp(-conj(s.exp))];
    if strcmpi(type, 'polres')
        r = @(th) scl*(eval_polres(th, pol, res, dom))+const;
        return
    else %type = 'zt'
        r = @(z) scl*(eval_zt(z, pol, res, dom))+const; 
        return
    end
    
%types where we need to sample and build an rfun:
elseif ( strcmpi(type, 'cpqr') || strcmpi(type, 'aaa') )
    dom = s.domain; a = dom(1); b = dom(2); 
    if isempty(nums)
        NN = s.res;
    else
        NN = nums; 
        NN = NN + mod(2,1); 
        NN = (NN -2)/2; %sample rate = 2*NN+2 
    end
    mmax = 2*length(s.exp);
    const = s.const; %normalize sample
    scl = s.scl; 
    s.const = 0; 
    s.scl = 1; 
    P = poles(s); 
  
    h = []; 
    if ~fft_flag
        h = ift(s, 'polres'); 
    end
    
    if isempty(pts) %autosample version
        happy = false; 
        while ~happy %we check for real-valued poles to determine happiness.
        [vals, pts] = sample_rfun(s, h, NN, dom);    
            if strcmpi(type, 'aaa')
                if mmax > 20 %this helps with being too strict about tol 
                %and thereby introducing spurious poles. 
                    r = rfun(scl*vals+const, pts, 'mmax', mmax, 'dom', [a, b], 'tol',tol);
                else
                    r = rfun(scl*vals+const, pts, 'dom', [a, b], 'tol', tol);
                end
                pol = r.poles; 
            else
            pp = 2; %default oversampling parameter for CPQR routine
            [pol, res, zj, fj, wj] = pronyaaa_poles(vals,pts,P,dom,pp);
            end  
        %check for poles on the interval:
        if all(abs(imag(pol)) > 1e-11) 
            happy = true; 
        elseif 2*NN+2 > 4*(s.res) % max sample size. Sampling isn't fixing the problem. 
            %loosen tol and as a last resort, we try standard AAA:
            r = rfun(scl*vals+const, pts, 'dom', [a, b], 'tol', min(1e-3, 1e2*tol));
            warning('efun:ift: Unable to construct a stable solution at specified tolerance. Tolerance reset to %e. Spurious poles may be present or the degree of the rational may change.', 1e2*tol);
            return
        end
        NNo = NN; 
        NN = 2*NN; %double sample and try again.
        end   
    else  % sample and call ift: 
        %h = ift(s, 'polres'); 
        vals = h(pts); 
        NNo = length(vals); 
        NNo = NNo + mod(2,1); 
        NNo = (NN -2)/2;
        if strcmpi(type, 'aaa')
            if mmax > 20 %this helps with being too strict about tol 
                %and thereby introducing spurious poles. 
                r = rfun(scl*vals+const, pts, 'mmax', mmax+2, 'dom', [a, b], 'tol',tol);
            else
                r = rfun(scl*vals+const, pts, 'dom', [a, b], 'tol', tol);
            end
        else %cpqr type: 
            pp = 2; % We try with O.S. parameter pp. It is increased to 4
            % inside pronyaaa_poles if spurious poles are detected. 
            [pol, res, zj, fj, wj] = pronyaaa_poles(vals,pts,P,dom,pp);
        end  
   end
   %assign r props here if needed:
   if strcmpi(type, 'cpqr')
      r = rfun(); 
      r.poles = pol; 
      r.residues = res; 
      r.nodes = zj; 
      r.vals = fj; 
      r.weights = wj; 
      r.tol = s.tol; 
      r.const = const; 
      r.scl = scl;
      r.res = 2*NNo+1; 
      r.domain = dom; 
   end
   
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END MAIN%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = eval_polres(theta, pol, res, dom)

%if domain is not [0, 1), adjust it: 
sh = size(theta); 
theta = theta(:); 
theta = (theta - dom(1))/(dom(2) - dom(1)); %map to [0, 1)

z = exp(-2*pi*1i*(theta)); z = z(:); 

C = bsxfun(@minus, z,pol.');  
C = 1./C; 
vals = C*res; 
%should always be real-valued: 
vals = real(vals); 
vals = reshape(vals, sh); 

end

function vals = eval_zt(z, pol, res, dom)
a = dom(1); 
b = dom(2); 
sh = size(z); 
z = z(:); 
z = flip(z); % we have coded for z = exp(-2*pi*1i*x), this adjusts for 
            %input where z = exp(2*pi*1i*x). 
%adjust z to account for domain change:
%z = z*exp(2*pi*1i*x); z = z.^(1/(b-a));
C = bsxfun(@minus, z,pol.');  
C = 1./C; 
vals = C*res;  
%if all input are on unit circle, kill tiny imag part: 
if all(abs(1-abs(z)) < 1e-14)
    vals = real(vals);
end
vals = reshape(vals, sh); 
end    
    
function [vals, pts] = sample_rfun(s, h, NN, dom)
% get a sample of r for AAA: 
if isempty(h)
    cutoff = NN+1; 
   while ( (cutoff+1 > NN) && cutoff < 9999 )
    NN = 2*NN;
    coeffs = feval(s,0:NN).'; m = max(abs(coeffs)); 
    cutoff = min(NN+1, standardChop(coeffs, m*1e-1*s.tol));% we want to make sure we accurately get sample
   end
   coeffs = coeffs(1:cutoff); 
   coeffs = [flip(conj(coeffs(2:end)));coeffs]; 
   coeffs = ifftshift(coeffs); 
   vals = length(coeffs)*ifft(coeffs); 
   pts = linspace(dom(1),dom(2), 2*cutoff); 
   pts = pts(1:end-1).';
   vals = real(vals); 
else  
   pts = linspace(dom(1), dom(2), 2*NN+2); pts = pts(1:end-1).';
   vals = h(pts); 
end

end


  




