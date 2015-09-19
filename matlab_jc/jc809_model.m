% model: computes the model function and its derivatives with
% respect to the parameters.
% Should be re-defined for each model.
% In this case: A sum of Gaussians each having a variable
% position, amplitude, and width (after example 15.5.16 NR, 'fgauss').
% Note that since we call the model function once for all the data
% points 'dyda' is a 2-dim array with dimensions [ndata x na] rather
% than a 1-dim array with 'na' elements. Similarly, 'y' is an array
% of size 'ndata' and not a scalar.
% Input: x,a
% Output: y,dyda
function [y,dyda]=jc809_model(x,a)
na = length(a);
 y=0; % set all elements to 0
 for i = 1:3:na % loop over Gaussians
     % some coefficients
     ai = a(i); % amplitude
     ai1 = a(i+1); % position
     ai2 = a(i+2); % width
     ar = (x-ai1)/ai2;
     ex = exp(-ar*ar);
     fac = 2*ai*ex*ar;
     % the function and its derivatives
     y =y+ai*ex; % sum of Gaussians
     dyda(:,i) = ex; % column i : derivative - amplitude
     dyda(:,i+1) = fac/ai2; % column i+1: derivative - position
     dyda(:,i+2) = fac*ar/ai2; % column i+2: derivative - width
 end