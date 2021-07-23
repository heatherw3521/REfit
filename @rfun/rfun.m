classdef rfun
% r = rfun(f, x, varargin)
%
% RFUN class for representing f with a trigonometric rational 
% barycentric interpolant given by the object r. 
%
% We assume the function being represented is real-valued, 1D, and
% periodic or compactly supported.
%
% INPUT:
% f may be a function handle, vector of values, or an efun object. 
% 
% x is a set of evaluation points. If x is not supplied, it is assumed that
% x consists of equally spaced points over the domain. The number of points is 
% adaptively determined if f is a function handle or efun. Otherwise, 
% length(x) = length(f). 
%
% r = rfun(f, x, 'dom', [a b], 'tol', tol). Sets tol and domain. Default 
% domain = [0, 1], default tol = 1e-9. 
%
% r = rfun(f, x, 'deg', k) constructs a trigonometric rational interpolant 
% of type (k-1, k). 
%
%% SEE (ADD LITERATURE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASS CONSTRUCTOR:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = false)
        function r = rfun(varargin)
            % The main rfun constructor
            % Return an empty RFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                r.domain = [];
                return
            end         
            %main constructor (calls aaa_trig)
            r = constructor(r, varargin{:});           
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        poles
        residues   
        nodes    %interpolating points/barycentric nodes
        vals     %value at interpolating points
        weights  %barycentric weights
        const    %additive constant 
        scl      % normalization parameter
        domain
        res  %resolution parameter (# of points in domain discretization)
        tol
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %methods (Access = public, Static = false)      
    
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %methods
        
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = true, Hidden = true)
    %pronyaaa for trigonometric interpolation    
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa_trig(f, varargin); 
    %pronyaaa with prescribed poles (used for inverse Fourier transform of efuns)
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa_trig_poles(f, varargin); 
    end
  
end
