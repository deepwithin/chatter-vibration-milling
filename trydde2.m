clc
clear
close all

global m F omega_nx omega_ny xi_x xi_y

Kt = 600; % radial cutting coefficient (from data book) in MPa
Kr = 0.07; % tangential cutting coefficient (from data book)
xi_x = 0.035; % damping ratio in x dirn
xi_y = 0.035; % damping ratio in y dirn
omega_nx = 600; % natural frequency in Hz (x compnt)
omega_ny = 660; % natural frequency in Hz (y compnt)
kx = 5600*10^3; % stiffness in N/m
Ky = 5600*10^3; % stiffness in N/m
m = 0.001;
F = 1000;

a = 3; % axial depth of cut, ADOC in mm

N = 2000; % rpm
Z = 4; % no. of teeth on cutter
T = 60/(N*Z); % tooth passing period

lags = [T T];

tf = 0.12;

sol = dde23(@ddefunc, lags, @yhist, [0 tf]);

t = linspace(0,tf,7);
y = deval(sol, t);

figure
subplot(1, 2, 1);
plot(t, y);
grid on;

% subplot(1, 2, 2);
% plot3(y(1,:), y(2,:), y(3,:));
% grid on;

function yp = ddefunc(t, y, YL)
    
%     global m F omega_ny xi_y
    
    yl1 = YL(:,1); % lag on y
%     yl2 = YL(:,2); % lag on x
    y1 = y(1)
    y2 = y(2)
    
    YL
    
    yp = [ y2;
           -y1 - 1 + yl1(1) - y1 ] 
end

function y = yhist(t)

    y = [0; 0];

end