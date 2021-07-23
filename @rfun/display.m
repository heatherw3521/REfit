function display(X)
% Display information about an RFUN.
% DISPLAY(X) outputs important information about the RFUN X to the command
% window.
%
% It is called automatically when the semicolon is not used at the end of a
% statement that results in a RFUN.
%
% See also RFUN/DISP.

disp(' ');
disp([inputname(1), ' =']);
disp(' ');
disp(X);

end


