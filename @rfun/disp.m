function disp(F)
% Display an rfun to the command line.
%%
% See also rfun/display.

% deal with empty case.  
if ( isempty(F) )
    fprintf('    empty rfun\n')
    fprintf('\n');
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
len = length(F);                          % # of barycentric nodes
                             

% Display the information: 
disp(' rfun object')
%fprintf('\n');
fprintf( '   domain            # of poles \n');
fprintf( '  [%g,%4.2g]            %6i\n'...
        ,dom(1), dom(2),len);
fprintf('\n');

end

