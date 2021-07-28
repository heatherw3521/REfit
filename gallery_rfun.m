function varargout = gallery_rfun(name)
% Gallery function for rfuns
% 
% 'abs': abs(x)
%
% 'rand': quasi-random selection of clusters of poles.
%
% 'wells': an analytic function consisting of several well-like regions.
%
% 'wild': wild function from Trefethen (cite). 
%
% 'wilder': wild with some added singularities. 
%%

if ( nargin == 0 )
    names = {'abs', 'rand', 'wells', 'wild', 'wilder'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)
    
    case 'abs' 
        fa = @(x) abs(x-.5);
        n = 6000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = rfun(fa(x), x, 'tol', 1e-9); 
                 
    case 'rand'
        [f, fa] = gallery_efun('rand'); 
        f = ift(f); 
        
    case 'wells' 
        fa = @(x) wells(x); 
        n = 3000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = rfun(fa(x), x, 'tol', 1e-8); 
        
    case 'wild' 
        fa = @(x) wild(x); 
        n = 4000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = rfun(fa(x), x, 'tol', 1e-5); 
        
    case 'wilder' 
        fa = @(x) wilder(x); 
        n = 4000; x = linspace(0, 1, 2*n+2); x = x(1:end-1).';
        f = rfun(fa(x), x, 'tol', 1e-6); 
    %  error if the input is unknown.
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
