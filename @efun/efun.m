classdef efun
% % s = efun(f) 
%
% EFUN class for representing sums of weigthed decaying complex exponentials. 
%
% These sums can be used to approximate the Fourier series of f (default)
% or f directly (include flag 'values').
% 
% INPUT: 
% f is a vector of values [x_0, x_1, ...x_{2N}], where x_j = j/2N+1.
% f can also be a function handle, chebfun, or rfun object.
%
% s = efun(f, x,'coeffs'). f is a consecutive set of Fourier coefficients,
%  'x' is a vector of modes, by default x = [-N...0...N]. 
%
% s = efun(f, 'domain', [a b], 'tol', tol). Sets tol and domain. 
% domain is time/value space. 
%
% default tol = 1e-10*max(f), default domain = [0, 1].
%
% s = efun(f, varargin, 'values'): constructs an exponential sum s
% that approximates f directly from equally spaced samples on the given domain. 
%
%% SEE (literature).

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% CLASS CONSTRUCTOR:
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = false)
        function s = efun(varargin)
            % The main constructor
            
            % Return an empty EXPFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                s.domain = [];
                return
            end         
            %main constructor 
            s = constructor(s,varargin{:});           
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        exp   % lam, where terms in the sum are w*exp(lam*k)
        weights % w, where terms in the sum are w*exp(lam*k)
        domain %domain of the function f (inverse Fourier transform of s
               % if space = Fourier.)
        space % is the exponential sum representing a function
              % in Fourier space or time domain?
        const %constant added on to exponential sum
        scl   % scale for exponential sum
        res %resolution limit originally used to construct approximation
            % res=L means that s was constructed with a least squares fit
            % to Fourier coeffs [-L, L].
    end
    
    properties (Access = private)
        sv % singular values of the Hankel matrix used in 
           % Prony's method. 
        tol % orginal tolerance parameter for constructing efun.
    end
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true, Hidden = true )
        [z, ss] = coneigen(cp,fs,tol, hp) %regularized prony's method (RPM)
        [z, ss] = coneigen_rand(f, l, tol) % RPM with randomized SVD
        [r, ss, w] = prony_compress(cp,m,tol) %RPM rectangular 
        [w,r,ss, n]= vals2efun(f,x, varargin) %work on vals
        [w,r,ss, n]= coeffs2efun(f, varargin) %work on coeffs
        cutoff = standardChop(coeffs, tol) %chop for fourier coeffs
        w = trig_Baryweights(pts, dom) %barycentric weights for trig poly
        [coeffs, fvals, s_grid] = get_sample(f, space, dom, tol) %auto-sampling
    end
   
end
