function [smoothedtest,smoothedactual]=resampleCI(FFvals)
% [smoothedtest,smoothedactual]=resampleCI(FFvals,5);
% Making SEM bars with resampling
% This is based on a Science paper from back in the day
% See ScienceStatReview.pdf in my stats folder on desktop of old laptop (Fig5)

nv=length(FFvals);
smootherwidth=20; % vary this
smoothedactual=runningmedian(FFvals,smootherwidth); % smooth - vary this function?
resamplesize=1000;
clear smoothedtest
for i=1:resamplesize
    randvec1=sort(ceil(rand(1,nv)*nv));
    smoothedtest(i,:)=runningmedian(FFvals(randvec1),smootherwidth); % smooth 
end

% figure;plot(smoothedactual,'b')
% hold on;plot(prctile(smoothedtest,95),'k')
% hold on;plot(prctile(smoothedtest,5),'k')
% hold on;plot(prctile(smoothedtest,50),'k')   