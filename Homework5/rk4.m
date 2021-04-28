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