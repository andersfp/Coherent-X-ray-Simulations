function str = timestr(x)
% Takes the current time and date when called and formats it as a single
% string output. The formatting depends on the input.

% Get the current time and date
c = clock;

% Convert time and date numbers to strings
Y = num2str(c(1),'%04.0f');
M = num2str(c(2),'%02.0f');
D = num2str(c(3),'%02.0f');
h = num2str(c(4),'%02.0f');
m = num2str(c(5),'%02.0f');
s = num2str(c(6),'%02.0f');

% Format the output
switch x
    case 1
        str = [Y '_' M '_' D];
    case -1
        str = [D '_' M '_' Y];
    case 2
        str = [h '_' m '_' s];
    case 3
        str = [Y '_' M '_' D '_' h '_' m];
    case 4
        str = [Y '_' M '_' D '_' h '_' m '_' s];
    otherwise
        warning('Not recognized formatting code.');
        str = '';
end

