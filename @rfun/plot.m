function  plot(f,varargin)
% Basic plot for rfun object: plots f on its domain. 
%
% See also: rfun/plotcoeffs, rfun/poleplot. 
%%

% Initialization:
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
[F, X] = sample(f); 
plot(X, F, out{:}, lineStyle{:}, pointStyle{:})

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
% THIS CODE IS BASED ON PLOT CODES FROM CHEBFUN DEVELOPERS. 
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
%
% PARSEPLOTSTYLE   Parse inputs to PLOT. Extract 'lineWidth', etc.
%   [L, P, J, D, OTHER] = PARSEPLOTSTYLE(VARARGIN) parses the inputs VARARGIN and
%   strips out inputs to the MATLAB/PLOT() that should only be included once.
%   For example, 'LineWidth' or 'MarkerSize'. Those options which correspond to
%   the Line part of the plot are returned in L, those corresponding to the
%   discrete points are returned in P, options for the 'jumpLine' appear in J,
%   and all other inputs are returned in OTHER as a cell array.
%

lineOpts = {'LineStyle', 'LineWidth'};
pointOpts = {'Marker', 'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'};

% Initialise:
lineStyle = {};
pointStyle = {};

k = 1; % Look at all remaining arguments.
while ( k < numel(varargin) )
    vk = varargin{k};

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
        varargin(k:k+1) = [];
    else
        k = k + 1;
    end

end

% Assign the remaining arguments to OUT:
out = varargin;

end
