function interactive5D(src,event,rgb)
% Make plot of the comparison of data and reconstruction interactive.

% Get relevant sizes
n = size(rgb,4);
n_omega = size(rgb,5);

% Get the current index
j = str2double(src.Name);

% Get the scroll direction
c = event.VerticalScrollCount;

% Update the index
j = j - c;

% Check bounds of the index
if j < 1 || j > n_omega
    return
end

% Update the plot
for i = 1:n
    src.Children(i).Children.CData = rgb(:,:,:,n - i + 1,j);
end

% Update the index
src.Name = num2str(j);

