%% LT, see bk24bk68 analysis - this code needs work.

clear all;
% use below to get percentiles each day
freq_vals = evtaf_freq('batch.catch.keep', [2000 3500], 'd', 128, 'obs0', 0); % after labeling, only get actually triggered syl

ptiles=prctile(freq_vals(:,2),[5 30 50 70 95])
cov=std(freq_vals(:,2))/mean(freq_vals(:,2))
figure; hist(freq_vals(:,2),20);

hit_rate=triglabel2('batch.catch.keep','d','c','',1,0,0,0)
hitrate=sum(hit_rate)
