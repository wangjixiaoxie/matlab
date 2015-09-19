function [ angle ] = db_angle_2d_vectors(matrix1, matrix2 )
%db_angle_2d_vectors Summary of this function goes here
%   Detailed explanation goes here

x1 = matrix1(1);
y1 = matrix1(2);

x2 = matrix2(1);
y2 = matrix2(2);

angle = mod(atan2(x1*y2-x2*y1,x1*x2+y1*y2),2*pi);


end

