function [ CI_lo, CI_high, p ] = db_diff_sample_mean_boot(population1,population2,samplesize,boot_number,alpha)
%db_diff_sample_mean_boot Finds the difference in mean for a sample of your population.
%Does that boot_number of times. If your CI does not include 0, then it is
%significantly different. (Set alpha to a decimal, like 0.05)

if nargin < 3
    samplesize = 30;
    boot_number = 10000;
    alpha = .05;
elseif nargin < 4
    boot_number = 10000;
    alpha = .05;
elseif nargin < 5
    alpha = .05;
end

converted_alpha = alpha*100;


parfor i = 1:boot_number
    p(i) = mean(randsample(population1,samplesize))-mean(randsample(population2,samplesize));
end

CI_lo = prctile(p,converted_alpha/2);
CI_high = prctile(p,100-converted_alpha/2);


end

