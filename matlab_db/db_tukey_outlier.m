function [ matrix outliers_high outliers_low ] = db_tukey_outlier(matrix, which_column)
%db_tukey_outlier Gets rid of outliers in a matrix using Tukey's method
%   Given a matrix 'A', specify which column you are using to look for
%   outliers. It then calculates the Tukey min and max for outliers (Y > Q3
%   + 1.5*IQR and Y < Q1 - 1.5*IQR). Then it removes those entire rows from
%   the matrix.

if ~exist('which_column', 'var')
    which_column = 1;
else
end

outlier_max = (median(matrix(:,which_column))+iqr(matrix(:,which_column))/2) + 1.5*iqr(matrix(:,which_column));

outlier_min = (median(matrix(:,which_column))-iqr(matrix(:,which_column))/2) - 1.5*iqr(matrix(:,which_column));

outliers_high = find(matrix(:,which_column) >= outlier_max);
outliers_low = find(matrix(:,which_column) <= outlier_min);

matrix = removerows(matrix, 'ind', matrix(:,which_column) >= outlier_max);
matrix = removerows(matrix, 'ind', matrix(:,which_column) <= outlier_min);


end

