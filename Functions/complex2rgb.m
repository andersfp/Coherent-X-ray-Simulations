function RGB = complex2rgb(A,cmap)
% Convert the complex 2D image A into an RGB image with the intensity given
% by the amplitude of A and color given by the phase of A. The amplitude of
% A is taken in the range [0 1].

% Separate the amplitude and the phase
P = angle(A);
A = abs(A);

% Cap the amplitude at 1
A(A > 1) = 1;

% Convert the phase to indices
n = size(cmap,1);
P = round((P + pi)/(2*pi)*(n - 1)) + 1;

% Generate each color
RGB = A(:).*cmap(P(:),:);

% Reshape the output
s = size(A);
RGB = reshape(RGB,[s 3]);

