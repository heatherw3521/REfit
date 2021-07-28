%main constructor for rfuns: 
function r = constructor(r, op,varargin) 
%
% the main constructor for rfuns 
%%

% pass efuns to efun ift function: 
if isa(op, 'efun')
    r = ift(op, varargin{:}); 
    return
end

% parse input
[dom, tol, type, deg, samples, locs, poles, cleanup] = parse_input(op, varargin{:});

%%
% call appropriate pronyAAA 
switch type  
    case 'function_handle'
        % adaptive sampling + pronyAAA
                    %call pronyAAA: 
            if isempty(deg)
                deg = 100; % deg = max # iterations. 
            end
            [~, pol, res, zj, fj, wj, err, N, const, scl] = pronyaaa_auto(op, dom, tol, deg, cleanup);
            tol = err;      
    case 'samples'
        if ~isempty(poles)
            %call pronyAAA_poles
        else
            %call pronyAAA: 
            if isempty(deg)
                deg = min(100, length(samples)); % deg = max # iterations. 
            end
            % we want a mean zero function: 
            const = mean(samples); 
            samples = samples - const; 
            scl = max(samples); 
            %what to do about zero function:
            % we assign arbitrary poles/nodes/weights.
            if abs(scl) <1e-13 
                r.domain = dom;
                r.const = const;
                r.poles = [dom(1)+1i;dom(1)-1i];
                r.weights = [1; 1];
                r.nodes = [diff(dom)/4; 3*diff(dom)/4];
                r.vals = [0;0];
                r.scl = 1; 
                r.const = 0; 
                r.tol = tol; 
                r.res = 2; 
                return
            end   
            samples = samples/scl;
            [~, pol, res, zj, fj, wj, ~, N] = pronyaaa(samples, locs, dom, tol, deg, cleanup);
        end
end
%% assign properties    
r.poles = pol; 
r.residues = res; 
r.nodes = zj; 
r.vals = fj; 
r.weights = wj; 
r.domain = dom; 
r.res = N; 
r.tol = tol; 
r.const = const;
r.scl = scl; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [dom, tol, type, deg, samples, locs, poles, cleanup] = parse_input(op, varargin)      
% determine domain: 
  dom = [0, 1]; %default domain
  if isa(op, 'chebfun') 
     dom = op.domain; 
  elseif ~isempty(varargin)
     domid = find(strcmp(varargin,'dom'));
     if ~isempty(domid)
       dom = varargin{domid+1};  
     end   
     domid = find(strcmp(varargin,'domain'));
     if ~isempty(domid)
       dom = varargin{domid+1};  
     end
  end
  
  % determine tolerance
  tol = 1e-9; 
  tolid = find(strcmp(varargin,'tol'));
  if ~isempty(tolid)
     tol = varargin{tolid+1};  
  end
  
  %determine if a degree is given:
  deg = []; 
  degid = find(strcmp(varargin,'deg'));
  if ~isempty(degid)
     deg = varargin{degid+1};  
  end
  
  degid = find(strcmp(varargin,'mmax'));
  if ~isempty(degid)
     deg = varargin{degid+1};  
  end
  
  %adjust tol if degree is given and no tol was specified: 
  if ~isempty(deg) && isempty(tolid) 
    tol = 1e-11; 
  end
  
  %determine if poles are given: 
  polid = find(strcmp(varargin,'poles'));
  if ~isempty(polid)
     poles = varargin{polid+1};  
  else
     poles = []; 
  end
  
  if any(strcmpi(varargin, 'cleanup_off'))
      cleanup = 0; 
  else
      cleanup = 1; 
  end
  
  %determine type and get samples if needed: 
  if isa(op, 'double') %vector of samples provided
      samples = op(:); 
      if isa(varargin{1}, 'double') %evaluation points?
          locs = varargin{1};
          if min(locs) < dom(1) || max(locs) > dom(2)
              error('rfun:rfun:constructor: sample locations outside specified domain')
          end
      else
          N = length(samples); 
          locs = (dom(2)-dom(1))*(0:N-1)/N + dom(1);
          locs = locs.';
          if any( strcmpi(varargin,'coeffs'))
          % convert Fourier samples and call aaa_trig
            samples = ifft([flip(conj(samples(2:end))); samples]);
            locs = (dom(2)-dom(1))*(0:N-1)/N + dom(1);
          end
      end
      type = 'samples'; 
  elseif isa(op, 'function_handle')
      if ~isempty(varargin) && isa(varargin{1}, 'double') 
          if length(varargin{1}) > 1 %just sample at eval pts; 
              locs = varargin{1};
              if min(locs) < dom(1) || max(locs) > dom(2)
                error('rfun:rfun:constructor: sample locations outside specified domain')
              end
              samples = op(locs);   
          else                       %sample size
              N = varargin{1}; 
              locs = (dom(2)-dom(1))*(0:N-1)/N + dom(1);
              locs = locs.';
              samples = op(locs);
          end
       type = 'samples'; 
      else
       %adaptive sampling: 
       type = 'function_handle'; 
       locs = []; 
       samples = []; 
      end
  else
    error('rfun:constructor:cannot parse input')
  end
  
end
       
              
              

             
             

             