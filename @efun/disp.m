function disp(F)
% Display an efun to the command line.
% 
%%
% See also efun/display.

% deal with empty case.  
if ( isempty(F) )
    fprintf('    empty efun\n')
    fprintf('\n');
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
len = length(F);                          % # of terms in sum
spc = F.space; 


% Display the information: 
if strcmpi(spc, 'Fourier')   
    disp(' efun object in Fourier space')
    %fprintf('\n');
    fprintf( '   domain (assoc. rational)           length of sum \n');
    fprintf( '        [%g,%4.2g]                       %6i\n'...
        ,dom(1), dom(2),len);
else
    disp('   efun object \n ')
    fprintf('       domain             length of sum       \n');
    fprintf('      [%g,%4.2g]             %6i\n', dom(1), dom(2),len);
end
fprintf('\n');
end

