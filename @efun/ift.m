function r = ift(s, varargin)
% The inverse Fourier transform of s.  
%
% For s in Fourier space, r = ift(s) constructs an rfun representing a
% barycentric trigonometric rational corresponding to the inverse Fourier 
% transform of s (See [2]). 
%
% ift(s, 'pr') returns a handle to the pole-residue form of r(x).
% ift(s, 'zt') returns a handle for the pole-residue form of 
% r(z) = r(x), where z = exp(2*pi*1i*x). 
%
% For rfuns:
% ift(s, 'cpqr')  uses CPQR pivoting to do AAA (default) (See [2]).
%
% ift(s, 'aaa') applies standard AAA w/ greedy pivoting (See [1]).
%
% If the flag 'fft' is included, the sample used to construct r consists of 
% sampling s and then applying the ifft. Otherwise, the pole-residue form of 
% r is used for sampling. 
%
% ift(s, sample, locs), ift(s, 'type', sample, locs) constructs r using sample, 
% a vector of samples of r. 
%
% See also rfun/ft. 
%%

if isempty(s)
    r = []; 
    return
end

%% 
% determine type of ift:
type = 'cpqr';

if ~isempty(varargin)
    if ~isnumeric(varargin{1})
    type = varargin{1};
    varargin{1} = {}; 
    end
end

% deal with pole residue forms: 
if (strcmpi(type, 'polres')  || strcmpi(type,'zt') )
    dom = s.domain;
    s.domain = [0, 1]; %shift to get correct poles:
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
        r = @(z) scl*(eval_polres(z, pol, res, dom))+const; 
        return
    end
elseif ( strcmpi(type, 'cpqr') || strcmpi(type, 'aaa') )
% for these, we first need to sample r: 
    [a, b] = s.domain; 
    NN = s.res;
    mmax = length(s.nodes);
%parse input:
    vals = []; pts = []; 
    while ~isempty(varargin)% see if the fft flag is there
        if isnumeric(varargin{1}) 
            if any(size(varargin{1}) > 1) %sample provided
                fft_flag = 0; 
                vals = varargin{1}; vals = vals(:); 
                pts = varargin{2}; pts = pts(:); 
                varargin{1} = {}; 
                varargin{2} = {}; 
            else % sampling rate = varargin{1}. 
                NN = varargin{1};
                varargin{1} = {}; 
            end
        elseif ischar(varargin{1}) %use fft sampling
            fft_flag = 1; 
            varargin{1} = {}; 
        else
            error('efun:ift:cannot parse input')
        end
    end
        
    %get sample if needed
    if isempty(vals) 
        if fft_flag
            coeffs = feval(s,0:NN);
            coeffs = [flip(conj(coeffs(2:end)));coeffs]; 

            coeffs = ifftshift(coeffs); 
            vals = length(coeffs)*ifft(coeffs); 
            pts = linspace(0,1, 2*NN+2); 
            pts = pts(1:end-1); pts  = pts(:);
            vals = real(vals); 
        else
            h = ift(s, 'polres'); 
            pts = linspace(a, b, 2*NN+2); pts = pts(1:end-1).';
            vals = h(pts); 
        end
    end
    %%
    % now call appropriate ift: 
    if strcmpi(type, 'aaa')
        if mmax > 20
            r = rfun(vals, pts, 'mmax', 2*mmax+2);
        else
            r = rfun(vals, pts);
        end
    else %use CPQR pivoting
        pol = poles(s); 
        % set up matrix (denominator of Barycentric form): 
        V = bsxfun(@(poles, pts) cot(2*pi/(dom(2)-dom(1))*...
        (poles - pts)./2), pol, pts.');
        V = [V; vals.']; 
        [~, ~, idx] = qr(V, 'vector'); 
        pp = 2; %oversampling parameter
        idx = idx(1:2*mmax+pp); 
        kpts = pts(idx);
        %add row of vals: 
        kvals = vals(idx); 
        [~, polr, resr, zerr, zjr, fjr, wjr] = AAA_poles(vals,pts...
        ,pol,kpts,kvals);
        %TO DO: add in a test to evaluate accuracy/check for real-valued
        %poles. 
        %
        % if real-valued poles appear or if the accuracy is low, 
        % try sampling more densely or using a higher oversampling 
        % parameter. 
        %
        %shift poles to interval
        polr = mod(real(polr), dom(2)) + 1i*imag(polr); 
        r = rfun(); 
        r.baryweights = wjr; 
        r.barynodes = zjr; 
        r.baryvals = fjr; 
        r.poles = polr;
        r.residues = resr;
        r.zeros = zerr;
        r.resolve_param = length(pts); 
    end
else
    error('efun:ift:cannot parse input')
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = eval_polres(theta, pol, res, dom)
%z = exp(-2*pi*1i*theta); z = z(:);

%if domain is not [0, 1), adjust it: 
sh = size(theta); 
theta = theta(:); 
theta = (theta - dom(1))/(dom(2) - dom(1)); %map to [0, 1)

z = exp(-2*pi*1i*theta); z = z(:); 

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
%adjust z to account for domain change:
z = z*exp(2*pi*1i*a); z = z.^(1/(b-a));
C = bsxfun(@minus, z,pol.');  
C = 1./C; 
vals = C*res;  
vals = reshape(vals, sh); 
end    
    
  




