close all

y_init = [0.994; 0.0; 0.0; -2.00158510637908252240537862224];
tspan = [0.0, 17.1];
[t, y] = ode45(@threebody, tspan, y_init);
figure
plot(y(:,1), y(:, 3))
xlabel('u_1')
ylabel('u_2')
title('Solution using ode45')


steps = 1000;
h = (tspan(2) - tspan(1)) / steps;
[t_rk, y_rk] = rk4(@threebody, tspan, y_init, h);
figure
plot(y_rk(1,:), y_rk(3,:))

function [t,y] = rk4(f,tspan,y0,h)
y0 = y0(:); % make sure y0 is a column vector
m = length(y0); % problem size
t = tspan(1):h:tspan(2); % output abscissae
N = length(t)-1; % number of steps
y = zeros(m,N+1);
y(:,1) = y0; % initialize
% Integrate
for i=1:N
% Calculate the four stages
K1 = feval(f, t(i),y(:,i) );
K2 = feval(f, t(i)+.5*h, y(:,i)+.5*h*K1);
K3 = feval(f, t(i)+.5*h, y(:,i)+.5*h*K2);
K4 = feval(f, t(i)+h, y(:,i)+h*K3 );
% Evaluate approximate solution at next step
y(:,i+1) = y(:,i) + h/6 *(K1+2*K2+2*K3+K4);
end
end

function [ dydt ] = threebody( t, y )
%THREEBODY Summary of this function goes here
%   Detailed explanation goes here
    mu = 0.012277471;
    mu_h = 1 - mu;
    D1 = ((y(1) + mu) ^ 2 + y(3) ^ 2) ^ 1.5;
    D2 = ((y(1) - mu_h) ^ 2 + y(3) ^ 2) ^ 1.5;
    dydt(1) = y(2);
    dydt(2) = y(1) + 2 * y(4) - mu_h * (y(1) + mu) / D1 - mu * (y(1) - mu_h) / D2 ;
    dydt(3) = y(4);
    dydt(4) = y(3) - 2 * y(2) - mu_h * y(3) / D1 - mu * y(3) / D2;
    dydt = dydt';
end