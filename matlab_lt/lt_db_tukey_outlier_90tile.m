function [ matrix outliers_high outliers_low ] = db_tukey_outlier(matrix, which_column, numIQR)
%% LT 1/5/15 - fhcanged to use difference between 10th and 90th percentiles (instead of IQR as in tukey)

%% LT 11/12/14 - now can specific how many IQRs to add on both sides
%%
%db_tukey_outlier Gets rid of outliers in a matrix using Tukey's method
%   Given a matrix 'A', specify which column you are using to look for
%   outliers. It then calculates the Tukey min and max for outliers (Y > Q3
%   + 1.5*IQR and Y < Q1 - 1.5*IQR). Then it removes those entire rows from
%   the matrix.

if isempty(numIQR)
    numIQR = 1.5;
else
end


if ~exist('which_column', 'var')
    which_column = 1;
else
end

outlier_max = (median(matrix(:,which_column))+lt_inter90tile(matrix(:,which_column))/2) + numIQR*lt_inter90tile(matrix(:,which_column));

outlier_min = (median(matrix(:,which_column))-lt_inter90tile(matrix(:,which_column))/2) - numIQR*lt_inter90tile(matrix(:,which_column));

outliers_high = find(matrix(:,which_column) >= outlier_max);
outliers_low = find(matrix(:,which_column) <= outlier_min);

matrix = removerows(matrix, 'ind', matrix(:,which_column) >= outlier_max);
matrix = removerows(matrix, 'ind', matrix(:,which_column) <= outlier_min);


end

