function display(X)
%DISPLAY   Display information about an EFUN.
%   DISPLAY(F) outputs important information about the EFUN F to the command
%   window.
%
%   It is called automatically when the semicolon is not used at the end of a
%   statement that results in a EFUN.
%
% See also DISP.

disp(' ');
disp([inputname(1), ' =']);
disp(' ');
disp(X);

end


