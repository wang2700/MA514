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

