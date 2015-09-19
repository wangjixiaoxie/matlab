function [ CI_lo, CI_high, p ] = db_sample_cohen_d_boot(population1,population2,samplesize,boot_number,alpha)
%db_ranksum_boot Runs ranksum for a sample of your population. Gives p value
%to make a confidence interval using bootstrp
%   Detailed explanation goes here


parfor i = 1:boot_number
    p(i) = db_cohens_d(randsample(population1,samplesize),randsample(population2,samplesize));
end

CI_lo = prctile(p,alpha/2);
CI_high = prctile(p,100-alpha/2);


end

