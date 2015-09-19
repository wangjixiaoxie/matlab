function [xmax,ymax] = pinterp(xs, ys)

% does parabolically interpolates the peak given three points

x = [xs.^2, xs, [1;1;1]];
xi = inv(x);
abc =xi*ys;

xmax = -abc(2)/(2*abc(1));
ymax = abc(1)*xmax^2+abc(2)*xmax+abc(3);
return;