function [ d ] = db_cohens_d(Sample1,Sample2)
%db_cohens_d Calculating Cohen's d (effect size) for two samples. Must be
%vectors. Inputs: Sample1, Sample2
%   d = mean(sample1) - mean(sample2) / standard deviation
%   standard deviation = sqrt((n_1-1)*var_1)+(n_2-1)*var_2) / n_1 + n_2)


mean1 = mean(Sample1);
mean2 = mean(Sample2);

var1 = var(Sample1);
var2 = var(Sample2);

n1 = length(Sample1);
n2 = length(Sample2);

s = sqrt(((n1-1)*var1 + (n2-1)*var2) ./ (n1+n2));

d = (mean1 - mean2) ./ s;




end

