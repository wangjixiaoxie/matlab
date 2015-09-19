function [ bhat_coef, bhat_dist ] = db_battacharyya( sample1, sample2, numberofbins, tukey_outlier )
%db_battacharyya Input: sample1, sample2, numberofbins, tukey_outlier
%(optional). Calculates the battacharyya coefficient and distance for two
%samples.
%   Battacharyya coef = sigma(i=1:n) sqrt(sigma(sample1)i sigma(sample2)i)
%   Battacharyya dist = 1 - sqrt(battacharyya coef)

%Runs tukey_outlier function if you want it to trim your samples
if ~exist('tukey_outlier', 'var')
else
    sample1 = db_tukey_outlier(sample1);
    sample2 = db_tukey_outlier(sample2);
end


%Calculates the domain of the two samples
min_value = min([min(sample1) min(sample2)]);
max_value = max([max(sample1) max(sample2)]);

domain = linspace(min_value,max_value,numberofbins);

%counts the number of values that fall in each bin of domain
count1 = histc(sample1, domain);
count2 = histc(sample2, domain);

%normalizes the count vectors so that they sum to 1
normalized_count1 = count1./(length(sample1) + length(sample2));
normalized_count2 = count2./(length(sample1) + length(sample2));

%calculates the bhattacharyya coefficient

bhat_coef = 0;

for i = 1:numberofbins
    bhat_coef = bhat_coef +...
        sqrt(normalized_count1(i) * normalized_count2(i));
end

%calculates the bhattacharyya distance
bhat_dist = sqrt(1 - bhat_coef);


end

