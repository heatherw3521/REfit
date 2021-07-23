%main constructor for efuns: 
function F = constructor(f, op, varargin) 
%
% the main constructor for exponential sums (efuns). 
%% 
F = f;
[F.space, F.domain, tol, type, samples, locs, deg, pronytype]...
    = parse_input(op, varargin{:}); 

%set some default options: 
chop_on = 1; 

%TO DO: deal with 'degree'/'length' option.
% For now, throw error. 
if ~isempty(deg)
    error('efun:constructor: the construction of fixed-length exponential sums is not yet available.')
end

switch type
%% case 1: parameters for a sum are input
    case 'parameters' %order = (weights, exponents, domain, space)
        F.weights = op(:); 
        exp= varargin{1};
        F.exp = exp(:); 
        F.sv = []; 
        F.res = [];
        F.const = 0; 
        F.scl = 1; 
        F.tol = tol; 
        return
%% case 2: efun
    case 'efun'
        if ~(deg == length(op))
            F = compress(op, deg); 
        elseif tol > op.tol
            F = compress(op, 'tol', tol); 
        else
            F = op; 
        end
        return
%% case 3: rfun input    
    case 'rfun'
        F = ft(op, 'tol', tol, deg); %call Fourier transform
        return
%% case 4:construct efun in Fourier space from coeffs
    case 'coeffs'
        %shift and scale: 
        const = samples(1); 
        samples(1) = 0; 
        scl = max(abs(samples)); 
        samples = samples/scl; 
        [w,r, L, ss]= coeffs2efun(samples, locs, chop_on, tol,pronytype);
        res = L;
%% case 5: construct efun in time domain from samples
    case 'vals'
        %% TO DO %%
        error('efun:constructor: constructing efuns on values not yet available.')
        %const = mean(samples); 
        %samples = samples - const; 
        %scl = max(abs(samples)); 
        %samples = samples/scl; 
        %[w,r, L, ss]= vals2exp(samples, locs, chop_on, deg, tol);
        %res = 2*L+1;
 end
%% BUILD THE EFUN OBJECT: 
F.weights = w; 
F.exp = r;
F.const = const;
F.scl = scl;
F.res = res;
F.sv = ss; 
F.tol = tol;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% END MAIN PROGRAM %%%%%%%%%%%%%%%%%%%
%% PARSE INPUT FUNCTION/GETS SAMPLES FOR CONSTRUCTION OF EFUN

function [space, dom, tol, type, samples, locs, deg, pronytype] = parse_input(op, varargin)

% set defaults
space = 'Fourier';
dom = [0, 1]; 
tol = 1e-10; 
deg = []; 
pronytype = []; 

%check for time or values flag:
if any(strcmpi(varargin, 'time')) || any(strcmpi(varargin, 'values')) ||...
        any(strcmpi(varargin, 'value'))
    space = 'time';
end

%check for domain flag:
[~, domid] = find(strcmp(varargin,'domain'));
if ~isempty(domid)
    dom = varargin{domid+1}; %for special types rfun, chebfun, domain is inherited in code below
end

%check for tol flag:
[~, tolid] = find(strcmp(varargin,'tol'));
if ~isempty(tolid)
    tol = varargin{tolid+1}; 
end

%check for degree or length flag:
[~, degid] = find(strcmp(varargin,'deg'));
if ~isempty(degid)
    deg = varargin{degid+1}; 
end

[~, degid] = find(strcmp(varargin,'degree'));
if ~isempty(degid)
    deg = varargin{degid+1}; 
end

[~, degid] = find(strcmp(varargin,'length'));
if ~isempty(degid)
    deg = varargin{degid+1}; 
end
[~, pid] = find(strcmp(varargin,'pronytype')); 
if ~isempty(degid)
    pronytype = varargin{pid+1}; 
end

%% parse types
type = []; 

% frst deal with parameters of exponential sums
if any(strcmp(varargin,'parameters')) || any(strcmp(varargin,'efun'))
    %build directly from exponents and weights: 
    samples = []; 
    locs = []; 
    deg = length(varargin{1}); 
    type = 'parameters';
    return
end
%%
if isa(op, 'rfun') %deal with rfuns
    dom = op.domain; %inherit properties from the rfun
    samples = []; 
    locs = []; 
    deg = length(op.poles);
    type = 'rfun';
    return
%%    
elseif isa(op, 'efun') %deal with efuns
    dom = op.domain;  
    samples = []; 
    locs = [];  
    if ~isempty(varargin) && isnumeric(varargin{1})
        deg = varargin{1}; % for input: efun(s, k), where k =< than length(s). 
    else
        deg = length(op); 
    end
    type = 'efun'; 
    return
%%    
elseif isa(op, 'double') %deal with vectors of samples: 
    if strcmpi(space, 'Fourier') %build efun in Fourier space
        if any(strcmpi(varargin, 'coeffs')) %samples are Fourier coeffs
            samples = op(:); 
            N = length(samples);                
            if length(varargin)==1 %assume coeffs are of the form [-N...0...N]
               if ~mod(N,2) 
                   error('efun:constructor:specify the Fourier modes with efun(samples, modes,...)')
               else 
                   locs = (-(N-1)/2:(N-1)/2).'; 
               end
            else    
                locs = varargin{1}; locs = locs(:); %coeff indices. 
            end
            [idx,~] = find(locs>=0);  % keep only non-neg coeffs
            locs = locs(idx); 
            samples = samples(idx); 
            N = length(samples);
            dd = diff(locs);
            if ~(abs(sum(dd)-N+1) == 0)
                error('efun: constructor: coefficient modes must be consecutive integers.')
            end
            %pad with zeros if needed: 
            if ~mod(N,2)
                samples = [samples;0];
                %N = N+1; 
                locs = [locs; locs(end)+1]; 
            end 
        else %samples are taken in time domain
            samples = op(:);                      
            N = length(op);
            % check whether sample locations are given: If they are, 
            % we need to throw an error when they aren't on the correct
            % grid.
            if ~isempty(varargin)
                plocs = varargin{1};
                if isa(plocs, 'double')
                   if length(plocs)==length(op) 
                       if ~(min(plocs) >= dom(1) && max(plocs) <= dom(2))
                           if ~any(mod(plocs, 1))
                             error('efun: constructor: sample is not on domain. Did you mean to include "coeffs" flag?')  
                           else
                             error('efun: constructor: sample is not on domain.')
                           end
                       end  
                       locst = (dom(2)-dom(1))*(0:N-1)/N + dom(1); 
                       if abs(norm(plocs(:)-locst(:))) > 5*1e-15
                           error(['efun:constructor: sample locations on domain [a, b] ',...
                           'should be of the form x_j = (j-1)*(b-a)/N + a for a length N sample. ',...
                           'If this is not possible,try using rfun instead.'])
                       end
                   end
                end
            end   
        [samples, locs] = samples2coeffs(samples); %translate to coeffs
        end
        type = 'coeffs'; 
    elseif strcmpi(space, 'time') %build efun in time space
        samples = op(:); 
        N = length(samples);
        if ~mod(N,2)
            error('efun:constructor:An odd number of samples is required')
        end
        %check for grid of location
        if isempty(varargin{1}) || ~isa(varargin{1}, 'double')
           %assume it is equally spaced points on domain
           locs = (dom(2)-dom(1))/N*(0:N-1) + dom(1); 
        elseif length(varargin{1})==N
            %check if the locations are at equally spaced points. 
            locs = varargin{1}; 
            locs = locs(:); 
            if ~(abs(sum(diff(locs)/locs(1))-N) == 0)
                error('efun:constructor: sample locations must consist of equally-spaced points.',...
                    'Try rfun for samples on non-equally spaced points.')
            end
        end
        type = 'time_samples';
    end
%%     
elseif isa(op,'function_handle') %deal with function handles
     if ~isempty(varargin) && isa(varargin{1}, 'double') %is there a sampling grid or a number of samples to take?
            locs = varargin{1}; locs = locs(:); 
            if length(locs) == 1 
                N = locs ; 
                N = N + mod(N+1,2);  %N must be odd
                N = N + 2*mod((N-1)/2, 2);  %odd number of positive coeffs
                locs = (dom(2)-dom(1))/N*(0:N-1) + dom(1); 
                samples = op(locs); samples = samples(:);            
            elseif norm(diff(locs)-ones(length(locs)-1,1)*(locs(2)-locs(1)))>1e-13 %throw error for non-equispaced
                error('efun:constructor: sample locations must consist of equally-spaced points.',...
                    'Try rfun for samples on non-equally spaced points.')
            else
             samples = op(locs); samples = samples(:); 
            end
            if strcmpi(space, 'Fourier')
                [samples, locs] = samples2coeffs(samples); %translate to coeffs
                type = 'coeffs';
            else
                type = 'time_samples';
            end
     else
        %adaptive sampling: 
        [cfs, samples, locs] = get_sample(op, space, dom, tol);
        if strcmpi(space, 'Fourier')
            samples = cfs ; 
            locs = (0:length(cfs)-1).';
            type = 'coeffs'; 
        else
            type = 'time_samples';
        end
     end
%%
elseif isa(op, 'chebfun') %deal with chebfuns
    %grab the coeffs
    dom = op.domain; %use domain from chebfun object. 
    N = length(op); 
    N = N + mod(N+1, 2);
    N = N + 2*mod((N-1)/2, 2);  %odd number of positive coeffs
    if strcmpi(space, 'Fourier') %construct efun on Fourier coeffs of op. 
        samples = trigcoeffs(op, N);
        samples = samples((N-1)/2-1:end); %keep pos. coeffs only
        N = length(samples); 
        locs = (0:(N-1)/2).'; %pos coeffs
    else %construct efun in value space (time domain)
        locs = (dom(2)-dom(1))/N*(0:N-1) + dom(1);
        samples = op(locs); 
    end
%%    
else %what could it be?
    error('efun:constructor:parsing error.')
end

end
%%
function [s, locs] = samples2coeffs(samples)
 N = length(samples);
 %this choice of N pads with correct amount of zeros if needed.
 N = N + mod(N+1,2);  %N must be odd
 N = N + 2*mod((N-1)/2, 2);  %odd number of positive coeffs
 s = 1/(length(samples))*fft(samples, N, 1); 
 %coeffs correspond to modes [0..(N-1)/2]
 s = s(1:(N-1)/2+1);
 locs = (0:(N-1)/2).';
end
 
            

    
  
        
    


    