
%This script will be used to calculate binary cross-correlation

first_trial=1;
last_trial=70;
triallength=32000*11;
maxlag=16000;
total_corr=0;
%First loop through each trial
for ii=first_trial:last_trial
    
    indicesa=find(spk1ind(spksin1clust,2)==ii);
    indicesb=find(spk1ind(spksin2clust,2)==ii);
   
    %now I just need two vectors of spiketimes
    spktmsa=spk1ind(spksin1clust(indicesa,1));
    spktmsb=spk1ind(spksin2clust(indicesb,1));

    spktrna=zeros(triallength,1);
    spktrnb=zeros(triallength,1);
    spktrna(spktmsa,1)=1;
    spktrnb(spktmsb,1)=1;
    
    corr_out=xcorr(spktrna,spktrnb, maxlag);
    corr_out=corr_out*length(spktmsb);

    total_corr=total_corr+corr_out;
end