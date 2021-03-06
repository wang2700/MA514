clear;
close all;

x = linspace(0,1,1000);
c = 3/2;
e1 = exp(x)-(1+x);
e2 = exp(x)-(1+c*x);

max(abs(e1))
max(abs(e2))

figure
plot(x, e1);
title('e1')
xlabel('x')
ylabel('e1')

figure
plot(x, e2);
title('e2')
xlabel('x')
ylabel('e2')

A = [1/3, 1/4; 1/4, 1/5];
b = [.5; exp(1)-7/3];
c = A^-1*b;
e2 = exp(x)-(1+c(1)*x+c(2)*x.^2); 

figure
plot(x, e2);
title('e2')
xlabel('x')
ylabel('e2')

max(abs(e2))