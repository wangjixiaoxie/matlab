function [ p ] = db_ttest2_boot(population1,population2,samplesize)
%db_ttest2_boot Runs ttest2 for a sample of your population. Gives p value
%to make a confidence interval using bootstrp
%   Detailed explanation goes here

[h,p] = ttest2(randsample(population1,samplesize),randsample(population2,samplesize));


end

