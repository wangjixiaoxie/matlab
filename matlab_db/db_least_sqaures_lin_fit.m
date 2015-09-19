function [ vec ] = db_least_squares_lin_fit( matrix )
%db_least_squares_lin_fit Returns a,b in y = ax + b for best fit line
%through a set of data points. Input a nx2 matrix (x,y).
%   Minimization function of D = Sigma([y_i - (a*x_i + b)]^2)

Sum_x2 = sum(matrix(:,1).^2);
Sum_x = sum(matrix(:,1));
Sum_y = sum(matrix(:,2));
Sum_xy = sum(matrix(:,1).*matrix(:,2));
n = size(matrix,1);

M = [Sum_x2, Sum_x; Sum_x, n];
P = [Sum_xy; Sum_y];

vec = mldivide(M,P);

end

