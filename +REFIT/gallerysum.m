function varargout = gallerysum(name)
% Gallery function for efuns
% 
% 'abs': abs(x)
%
% 'chirp': a chirplet.
%
% 'rand': a rational function with clusters of poles that are selected quasi-randomly.
%
% 'spline': a cubic spline on [0, 1].
%
% 'sqrt': sqrt(sin(pi*x)).
%
% 'wells': an analytic function consisting of several well-like regions.
%
% 'wild': the WILD function from Computing Numerically with Functions, 
% Trefethen 2007, rescaled to [0, 1]. 
%
% 'wilder': a function similar to WILD, but with singularities.
%
%%

if ( nargin == 0 )
    names = {'wild', 'spline', 'abs', 'chirp', 'sqrt', 'wilder', 'wells',...
        'rand'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)

    case 'wild' %wild from Trefethen
        fa = @(x) wild(x); 
        n = 6000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x), 'tol', 1e-11); 
    
    case 'wilder' %wild from Trefethen with added singularities
        fa = @(x) wilder(x); 
        n = 8000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x), 'tol', 1e-11); 
        
    case 'spline' %cubic spline
        f1 = @(x) 3/4*abs(x).^3 - 3/2*(x).^2 + 1; 
        f2 = @(x) 1/4*(2 - abs(x)).^3; 
        f3 = @(x) f1(x).*(abs(x) <= 1) + f2(x).*(abs(x) > 1 & abs(x) <= 2); 
        fa = @(x) f3(6*x-3)-.25; 
        n = 500; x = linspace(0, 1, 2*n+2); x = x(1:end-1).'; 
        f = efun(fa(x), 'tol', 1e-11); 
        
    case 'abs' 
        fa = @(x) abs(x-.5);
        n = 6000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x)-mean(fa(x)), 'tol', 1e-9); 
        
    case 'chirp'
        f = @(x) sin(73*pi*x.^2) + sin(29*pi*x);
        w = @(x) erfc(50*(x-0.8)).*erfc(-50*(x-0.2))/4;
        fa = @(x) (f(x)).*(w(x));%-0.003473067797593;
        n = 1000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x), 'tol', 1e-9); 
        
    case 'sqrt'
        fa = @(x) sqrt(sin(pi*x)); 
        n = 6000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x)-mean(fa(x)), 'tol', 1e-9);
        
    case 'rand'
        n = 10; rn = 10; clear poles
        poles = .5*((1:n).').^(-sqrt(2)); poles = .5 +2.^(-(n+1:2*n).') - (poles*1i);
        poles = [poles; rand(rn, 1) - 1i*rand(rn, 1)]; 
        nodes = -poles*2*pi*1i; 
        % .5 - 1i*poless];
        weights = (-1).^randi([1,2],rn+n,1).*(rand(rn+n,1) + rand(rn+n,1)*1i);
        f = efun(weights, nodes, [-1000 1000], 'parameters'); 
        f.space = 'Fourier'; 
        f.tol = 1e-8; 
        fa = @(x) feval(x, ss, 'values'); 
        
    case 'wells' 
        fa = @(x) wells(x); 
        n = 3000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = efun(fa(x), 'tol', 1e-6); 
    % Raise an error if the input is unknown.
    otherwise
        error('EFUN:GALLERY:unknown:unknownFunction', ...
            'Unknown function.')
end

% Only return something if there is an output argument.
if (nargout ==2)
    varargout{1} = f;
    varargout{2} = fa; 
elseif ( nargout ==1)
    varargout{1} = f; 
else
    % Otherwise, plot the function.
    plot(f)
    title([name ', length = ' num2str(length(f))])
    if ( ~isempty(ylim) )
        ylim(ylim)
    end
    shg
end

end




function s = wild(x)
% The 'wild' function from Computing Numerically with Functions, Trefethen 2007.
% rescaled to interval [0, 1). 

sz = size(x); 
x = x(:); 
f =  sin(pi*(2*x-1));
s = f;
for j = 1:15
    f = (3/4)*(1 - 2*f.^4);
    s = s + f;
end 
s = reshape(s, sz); 
end

function s = wilder(x)
% The 'wild' function from Computing Numerically with Functions, Trefethen 2007.
% rescaled to interval [0, 1) and with added singularities
 
sz = size(x); 
x = x(:); 
f =  sin(pi*((2*abs(x-.5)-1)));
s = f;
for j = 1:15
    f = (3/4)*(1 - 2*f.^4);
    s = s + f;
end 
s = reshape(s, sz); 
end

function f = wells(x)
% a function with well-like regions. 
sz = size(x); 
x = x(:); 
f =  sin(pi*(2*x-1));
for j = 1:6
    f = (3/4)*(1 - 2*f.^4);
end
f = f - mean(f); 
f = reshape(f, sz); 
end
