
function g = gaussian(x,mn,stdev)

g = (1/(stdev*sqrt(2*pi)))*exp(-0.5*((x-mn)/stdev).^2);
g=g/max(g);
