function display(X)
% Display information about an efun.
%   display(X) outputs important information about the efun X to the command
%   window.
%
%   It is called automatically when the semicolon is not used at the end of a
%   statement that results in an efun.
%
%%
% See also efun/disp.

disp(' ');
disp([inputname(1), ' =']);
disp(' ');
disp(X);

end


