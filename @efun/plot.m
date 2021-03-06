function  plot(f,varargin)
% Basic  plot for efun object:plots the associated trigonometric 
% rational ift(s) on itsdomain. 
%
%%
% See also: efun/plotcoeffs. 


%% Initialization:

if isempty(f)
    return %nothing to do. 
end

% Store the hold state of the current axis:
holdState = ishold;


% Remove global plotting options from input arguments.

if ~isempty(varargin)
[lineStyle, pointStyle, out] = ...
    parsePlotStyle(varargin{:});
else %assign a default style
    lineStyle = {'linewidth', 1.5};
    pointStyle = {}; 
    out = {}; 
end

%% create plot:
space = f.space; 
a = f.domain(1); b = f.domain(2); N = f.res; 
if strcmpi(space , 'Fourier')
    h = ift(f, 'polres'); 
    pts = linspace(a, b, 2*N+2); pts = pts(1:end-1).';
    plot(pts, h(pts), out{:},lineStyle{:}, pointStyle{:} )
else
    pts =linspace(a,b, N); 
    plot(pts, feval(f, pts), out{:}, lineStyle{:}, pointStyle{:} )
end


%% Misc:

% Return hold state to what it was before
% and set some basic style parameters:

if ( ~holdState )
    hold off
    set(gcf,'color','w')
    set(gca, 'fontsize', 18)
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%
%%
function [lineStyle, pointStyle,  out] = parsePlotStyle(varargin)
%PARSEPLOTSTYLE   Parse inputs to PLOT. Extract 'lineWidth', etc.
%   [L, P, J, D, OTHER] = PARSEPLOTSTYLE(VARARGIN) parses the inputs VARARGIN and
%   strips out inputs to the MATLAB/PLOT() that should only be in cluded once.
%   For example, 'LineWidth' or 'MarkerSize'. Those options which correspond to
%   the Line part of the plot are returned in L, those corresponding to the
%   discrete points are returned in P, options for the 'jumpLine' appear in J,
%   and all other inputs are returned in OTHER as a cell array.
%
% THIS CODE IS BASED ON CODE FROM CHEBFUN DEVELOPERS. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

lineOpts = {'LineStyle', 'LineWidth'};
pointOpts = {'Marker', 'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'};


% Initialise:
lineStyle = {};
pointStyle = {};

k = 1; % Look at all remaining arguments.
while ( k < numel(varargin) )
    vk = varargin{k};

    % Using strncmp() here is technically incorrect, but it's easier than
    % writing the code to properly do the name matching, and it's highly
    % unlikely that it will ever actually cause us any problems.
    if ( any(strncmpi(vk, lineOpts, 5)) )
        % Line option:
        lineStyle = [lineStyle, vk, varargin{k+1}]; 
        varargin(k:k+1) = [];
        
    elseif ( any(strncmpi(vk, pointOpts, 6)) )
        % Point option:
        pointStyle = [pointStyle, vk, varargin{k+1}];
        varargin(k:k+1) = [];
        
    elseif ( strcmpi(vk, 'color') )
        % Option for all:
        lineStyle = [lineStyle, vk, varargin{k+1}];
        pointStyle = [pointStyle, vk, varargin{k+1}];
        jumpStyle = [jumpStyle, vk, varargin{k+1}];
        varargin(k:k+1) = [];
    else
        k = k + 1;
    end

end

% Assign the remaining arguments to OUT:
out = varargin;

end


function [styleArgs, isValidStyle] = parseMATLABLineStyle(str)
%PARSEMATLABLINESTYLE   Validate and break a MATLAB line style into its pieces.
%   STYLEARGS = PARSEMATLABLINESTYLE(STR) takes the MATLAB line style
%   specification in STR and converts it into a cell array STYLEARGS of
%   'LineStyle', 'Marker', and 'Color' keyword pairs suitable for passing
%   to built-in PLOT.
%
%   [STYLEARGS, ISVALIDLINESTYLE] = PARSEMATLABLINESTYLE(STR) additionally
%   returns a logical variable that is TRUE when the style supplied in STR is a
%   well-formed MATLAB line style and FALSE otherwise.

    styleArgs = {};

    % Line style.  (This MUST be done first for '-.' to be parsed correctly.)
    [indS, indE, ~, line] = regexp(str, '--|-\.|-|:');
    str(indS:indE) = [];
    if ( ~isempty(line) )
        styleArgs = [styleArgs, 'LineStyle', line];
    end

    % Marker style.
    [indS, indE, ~, marker] = regexp(str, '[.ox+*sdv^<>ph]');
    str(indS:indE) = [];
    if ( ~isempty(marker) )
        styleArgs = [styleArgs, 'Marker', marker];
    end

    % Point style.
    [indS, indE, ~, color] = regexp(str, '[bgrcmykw]');
    str(indS:indE) = [];
    if ( ~isempty(color) )
        styleArgs = [styleArgs, 'Color', color];
    end

    isValidStyle = isempty(str);
end


