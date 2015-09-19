% lt modified so that randperm only takes one argument. m (my version of
% matlab doesn't allow two arguments?)

function [ CI_median, CI_lo, CI_high, sample_boot ] = lt_db_sample_boot_stat(population1, population2, func, samplesize, numofboot, alpha )
%db_sample_boot_stat Inputs: population1, population2, function,
%samplesize, numberofbootstrap, alpha. function must have the form
%@function (i.e. @mean)
%   Takes two large populations (population1 & population2) and performs 
%   the test function (@func) of interest for a random sample (N = samplesize)
%   numofboot times, and gives CI set by alpha.
%
%   example input(pop1,pop2,@(x,y)mean(x)-mean(y),30,10000,0.05) will calculate the
%   difference in the sample means with a sample size of 30 taken randomly
%   from pop1 and pop2, and do this 10,000 times.

sample_boot = zeros(numofboot,1);

for i = 1:numofboot
    randperm_population1_indices=randperm(length(population1));
    randperm_population2_indices=randperm(length(population2));
    sample_boot(i) = feval(func, population1(randperm_population1_indices(1:samplesize)),...
        population2(randperm_population2_indices(1:samplesize)));
end

CI_lo = prctile(sample_boot,100*(alpha/2));
CI_high = prctile(sample_boot,100*(1-alpha/2));
CI_median = median(sample_boot);

end

