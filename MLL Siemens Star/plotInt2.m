function plotInt2(x1,y1,x2,y2)
% Plot one column of y at a time, using the keyboard keys to browse through
% the columns.

% Set parameters
i = 1;
n = size(y1,2);

% Make the initial plot
plot(x1,y1(:,i),x2,y2(:,i));
set(gcf,'KeyPressFcn',@callback)
title(['i = ' num2str(i)]);

% Make the callback function
function callback(src,event)
    % Change the index on key press
    if strcmp(event.Key,'rightarrow')
        if i < n
            i = i + 1;
        end
    elseif strcmp(event.Key,'leftarrow')
        if i > 1
            i = i - 1;
        end
    end
    
    % Update the plot data
    src.Children.Children(1).YData = y2(:,i);
    src.Children.Children(2).YData = y1(:,i);
    title(['i = ' num2str(i)]);
    
end

end
