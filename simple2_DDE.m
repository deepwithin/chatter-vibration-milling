clc
clear
close all

lags = [1 0.5];

tf = 5;

sol = dde23(@ddefunc, lags, @yhist, [0 tf]);

t = linspace(0, tf, 100);
y = deval(sol, t);

figure
subplot(1, 2, 1);
plot(t, y);
grid on;

subplot(1, 2, 2);
plot3(y(1,:), y(2,:), y(3,:));
grid on;

function yp = ddefunc(~, y, YL)
    
    yl1 = YL(:,1); % lag: 1
    yl2 = YL(:,2); % lag: 0.5
    YL
    y
    yp = [yl1(1)
          y(1) - yl1(1) + yl2(2)
          y(2) - y(3)]

end

function y = yhist(t)

    y = [1 0 -1]';

end