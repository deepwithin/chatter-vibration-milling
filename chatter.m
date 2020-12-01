%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%_CHATTER VIBRATION ANALYSIS IN THIN WALL MILLING_%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
% format
format compact

%% Variable declaration

global omega_nx omega_ny xi_x xi_y Kt Kr kx Ky a N Z T theta time theta_entry theta_exit mx my fz f
global f1 c1 k1 f2 c2 k2 f1_sum f2_sum state

Kt = 600; % radial cutting coefficient (from data book) in MPa
Kr = 0.07; % tangential cutting coefficient (from data book)
xi_x = 0.035; % damping ratio in x dirn
xi_y = 0.035; % damping ratio in y dirn
omega_nx = 600; % natural frequency in Hz (x compnt)
omega_ny = 660; % natural frequency in Hz (y compnt)
kx = 5600*10^3; % stiffness in N/m
Ky = 5600*10^3; % stiffness in N/m

mx = 10; % modal mass
my = 10;

a = 3; % axial depth of cut, ADOC in mm

N = 2000; % rpm
Z = 4; % no. of teeth on cutter
T = 60/(N*Z); % tooth passing period
f = 1; % feed
fz = f/Z; % feed per tooth

lags = [T T]; % time lags for DDE equal to tooth passing period

tf = 3; % total time of simulation

theta = 0; % instantaneous theta initialized
theta_entry = 0 *pi/180; % in rads
theta_exit = 90 *pi/180; % in rads
x=0;
y=0;

state = 0; % represents time in steps of tooth passing time

time = linspace(0,tf,10000); % in secs

%% First calculation

f1_sum = 0;
f2_sum = 0;
fprintf('################ Simulation begins ################');
for tooth = 1:Z
        
    theta = 2*pi*N*state/60 + Z*2*pi/N;
    fprintf('< %d >-----',tooth);
    engagement = g(theta);
    f1 = ( Kt*Kr*a*(-sin(theta)) + Kt*a*(-cos(theta)) )*engagement;
    f2 = ( Kt*Kr*a*(-cos(theta)) + Kt*a*(sin(theta)) )*engagement;

    f1_sum = f1 + f1_sum;
    f2_sum = f2 + f2_sum;

end

c1 = -2*xi_x*omega_nx;
k1 = -omega_nx^2;

c2 = -2*xi_y*omega_ny;
k2 = -omega_ny^2;


options = ddeset('Event',@ChtrEvents); % custom trigger event is set
state = state + T
sol = dde23(@ddefunc, lags, @yhist, [0 tf], options);

%% Simulation calculations till final time

while sol.x(end) < tf
    state = state + T;
    fprintf('\n______________Integration Restart at %5.6f_____________\n', sol.x(end));
    fprintf('state value after this cycle = %f \n',state);
    
    f1_sum = 0;
    f2_sum = 0;
    
    for tooth = 1:Z
                
        theta = 2*pi*N*state/60 + (tooth-1)*2*pi/Z;
        fprintf('< %d >-----',tooth);
        engagement = g(theta);
        f1 = ( Kt*Kr*a*(-sin(theta)) + Kt*a*(-cos(theta)) )*engagement;
        f2 = ( Kt*Kr*a*(-cos(theta)) + Kt*a*(sin(theta)) )*engagement;
        
        f1_sum = f1 + f1_sum;
        f2_sum = f2 + f2_sum;
        
    end
    
    c1 = -2*xi_x*omega_nx;
    k1 = -omega_nx^2;
    
    c2 = -2*xi_y*omega_ny;
    k2 = -omega_ny^2;
    
    % calculating the solution
    sol = dde23(@ddefunc, lags, sol, [sol.x(end) tf], options);
    %y = deval(sol, t); % t is diif from time
    
end
fprintf('\n******************Simulation Done******************\n\n');

%% Roughness calculations

Ray = 0; % avg. roughness y dirn
for i = 1:length(sol.x)
    itr = i-1;
    if itr == 0
        prevTime = 0;
    else
        prevTime = sol.x(1,itr);
    end
    Ray = Ray + abs(sol.y(1,i))*( sol.x(1,i) - prevTime )/tf ;
end

fprintf('Average Roughness in Y direction = %f mm\n',Ray*1000);

Rax = 0; % avg. roughness x dirn
for i = 1:length(sol.x)
    itr = i-1;
    if itr == 0
        prevTime = 0;
    else
        prevTime = sol.x(1,itr);
    end
    Rax = Rax + abs(sol.y(3,i))*( sol.x(1,i) - prevTime )/tf ;
end

fprintf('Average Roughness in X direction = %f mm\n\n',Rax*1000);

Rqy = 0; % rms roughness y dirn
for i = 1:length(sol.x)
    itr = i-1;
    if itr == 0
        prevTime = 0;
    else
        prevTime = sol.x(1,itr);
    end
    Rqy = Rqy + ( (sol.y(1,i))^2 )*( sol.x(1,i) - prevTime )/tf ;
end
Rqy = sqrt(Rqy);
fprintf('Root mean square Roughness in Y direction = %f mm\n',Rqy*1000);

Rqx = 0; % rms roughness x dirn
for i = 1:length(sol.x)
    itr = i-1;
    if itr == 0
        prevTime = 0;
    else
        prevTime = sol.x(1,itr);
    end
    Rqx = Rqx + ( (sol.y(3,i))^2 )*( sol.x(1,i) - prevTime )/tf ;
end
Rqx = sqrt(Rqx);
fprintf('Root mean square Roughness in X direction = %f mm\n\n',Rqx*1000);

Rty = max(sol.y(1,:)) - min(sol.y(1,:));
Rtx = max(sol.y(3,:)) - min(sol.y(3,:));

fprintf('Total height of profile in Y direction = %f mm\n',Rty*1000);
fprintf('Total height of profile in X direction = %f mm\n\n',Rtx*1000);

%% Amplitude plots- Separate graphs

% figure
% plot(sol.x,sol.y(1,:));
% title('Vibration profile');
% % set(gca,'FontSize', 12);
% ylabel('y displacement in mm');
% xlabel('time');
% 
% figure
% plot(sol.x,sol.y(3,:));
% title('Vibration profile');
% % set(gca,'FontSize', 12);
% ylabel('x displacement in mm');
% xlabel('time');

%% Amplitude plots- Combined graph

figure
subplot(2,1,1)
plot(sol.x,sol.y(1,:));
title('Amplitude variation in Y direction')
set(gca,'FontSize', 14);
ylabel('y displacement in m');
xlabel('time');

subplot(2,1,2)
% plot(t, x, 'red');
plot(sol.x,sol.y(3,:));
title('Amplitude variation in X direction')
set(gca,'FontSize', 14);
ylabel('x displacement in m');
xlabel('time');

%% Function declarations
function yp = ddefunc(t, y, YL)
    
     global theta c1 k1 f1_sum fz mx c2 k2 f2_sum my
    
    yl1 = YL(:,1); % lag on y
    yl2 = YL(:,2); % lag on x
    y1 = y(1);
    y2 = y(2);
    x1 = y(3);
    x2 = y(4);

    yp = [ y2;
           c2*y2 + k2*y1 + f2_sum*fz*sin(theta)/my + f2_sum*sin(theta)/my*(x1 - yl2(2)) + f2_sum*cos(theta)/my*(y1 - yl1(1)) ;
           x2;
           c1*x2 + k1*x1 + f1_sum*fz*sin(theta)/mx + f1_sum*sin(theta)/mx*(x1 - yl2(2)) + f1_sum*cos(theta)/mx*(y1 - yl1(1)) ] ; 
end

function y = yhist(t)
    y = [0 0 0 0]';
end

function [value,isterminal,direction] = ChtrEvents(t,y,YL)
    global state
%     t
    if state-t == 0
        fprintf('Event triggered. Integration terminated\n');
    end
    value = [state - t; state - t; state - t; state - t] ;
    isterminal = [1; 1; 1; 1];
    direction = [0; 0; 0; 0];
    
end

function intermittent_check = g(theta)
    
    global theta_entry theta_exit
    fprintf('theta = %5.10f\n',theta);
    fprintf('theta in 0 to 360: %4.4f \n', mod(theta,2*pi));
    
    if (theta_entry < mod(theta,2*pi)) && (mod(theta,2*pi) < theta_exit)
        intermittent_check = 1;
        fprintf('tooth engaged\n');
    else
        intermittent_check = 0;
        fprintf('tooth is free\n');
    end
end

