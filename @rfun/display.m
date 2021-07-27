function display(X)
% Display information about an rfun.
% display(X) outputs important information about the rfun X to the command
% window.
%
% It is called automatically when the semicolon is not used at the end of a
% statement that results in a rfun.
%
% See also rfun/disp.

disp(' ');
disp([inputname(1), ' =']);
disp(' ');
disp(X);

end


