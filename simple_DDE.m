clc
clear
close all

lags = 1;

tf = 5;

sol = dde23(@ddefunc, lags, @yhist, [0 tf]);

t = linspace(0, tf, 100);
y = deval(sol, t);

figure
plot(t,y)

function yp = ddefunc(t, y, YL)

    y
    YL
    yp = 2*y - YL - 1

end

function y = yhist(t)

    y = 1;

end