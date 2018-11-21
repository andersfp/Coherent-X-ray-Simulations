% Initialization
clear;
close all;
clc;


%% Load the data
% File name
f = 'Falch_Pink_Resolution.txt';

% Load file
dat = dlmread(f,'\t');

% Extract the data
[x,ii,~] = unique(dat(:,1));
y = dat(ii,2);

% Plot the data
figure;
plot(x,y);


%% Separate the data
% Set each segment
x1 = [0.1 4.8 10.8 17.6 23.5 30.1 35.5 41.2 47.3];
x2 = [7.2 14.8 19.5 27.7 31.7 39.2 44.5 51.5 56.2];

% Find the indices
[~,i1] = min(abs(x - x1));
[~,i2] = min(abs(x - x2));

% % Get the top and bottom data
% a = [];
% b = [];
% for i = 1:4
%     a = [a;y(i1(2*i+1):i2(2*i))];
%     b = [b;y(i1(2*i):i2(2*i-1))];
% end
% a = mean(a);
% b = mean(b);
% 
% % Plot the data
% figure;
% plot(x,y,[x(1) x(end)],[a a],[x(1) x(end)],[b b]);

% Fit data
fun = @(a,b,c,sig,x) a.*erf((x - b)./(sqrt(2).*sig)) + c;
n = length(x1);
a = zeros(n,1);
ae = a;
b = a;
be = a;
c = a;
ce = a;
sig = a;
sige = a;
for i = 1:n
    ff = fit(x(i1(i):i2(i)),y(i1(i):i2(i)),fun,'StartPoint',[(-1)^(i)*0.07 mean(x([i1(i) i2(i)])) 0.92 0.4]);
    figure;
    plot(ff,x(i1(i):i2(i)),y(i1(i):i2(i)));
    fe = diff(confint(ff));
    a(i) = ff.a;
    b(i) = ff.b;
    c(i) = ff.c;
    sig(i) = 1e3*ff.sig;
    ae(i) = fe(1);
    be(i) = fe(2);
    ce(i) = fe(3);
    sige(i) = 1e3*fe(4);
end

% Plot the results
figure;
errorbar(1:n,abs(a),ae);
figure;
errorbar(1:n,b,be);
figure;
errorbar(1:n,c,ce);
figure;
errorbar(1:n,sig,sige);


