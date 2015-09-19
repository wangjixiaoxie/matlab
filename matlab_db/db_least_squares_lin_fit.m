function [ vec ] = db_least_squares_lin_fit( x,y )
%db_least_squares_lin_fit Returns a,b in y = ax + b for best fit line
%through a set of data points. Input a nx2 matrix (x,y).
%   Minimization function of D = Sigma([y_i - (a*x_i + b)]^2)

Sum_x2 = sum(x.^2);
Sum_x = sum(x);
Sum_y = sum(y);
Sum_xy = sum(x.*y);
n = size(x,1);

M = [Sum_x2, Sum_x; Sum_x, n];
P = [Sum_xy; Sum_y];

vec = mldivide(M,P);

end

