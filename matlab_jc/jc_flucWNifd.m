function [autocorrelations,widths,xpsd,psds]=jc_flucWNifd
% for 50 white noise files
% generate pitch curves within a window of width 200Hz to 800Hz
% NOTE: it doesn't appear that windowwidth or windowmean make a difference
% vary sigma and measure the zero-crossing point of the mean-subtracted
% autocorrelation figures.

    fs = 32000; % cbin sampling rate
    sigma=[0.5 1 2];
    
for kk=1:3
% Do this for 50 random wn files, 250ms in length
for jj=1:50
    wndata=rand(8000,1)*2-1;
    pitch=jc_PitchData208(wndata',1024,1020,sigma(kk),2600,3400,1,'obs0');
    pitches(jj,:)=pitch;
    autocorr=xcorr(pitch-mean(pitch));
    autocorrstor(jj,:)=autocorr;
    %%%%% This code finds the width between zero crossings, a metric of the
    %%%%% maximal width of artificially induced correlations by the method.
        for i=1:length(autocorr)-1
            increasing(i)=autocorr(i)<0&autocorr(i+1)>0;
            decreasing(i)=autocorr(i)>0&autocorr(i+1)<0;
        end
        goingup=find(increasing==1);
        goingdown=find(decreasing==1);
        [peak,peakindex]=max(autocorr);
        [p,rightind]=min(abs(goingdown-peakindex));
        [q,leftind]=min(abs(peakindex-goingup));
        leftpt=goingup(leftind)+abs(autocorr(goingup(leftind)))/(abs(autocorr(goingup(leftind)))+abs(autocorr(goingup(leftind)+1)));
        rightpt=goingdown(rightind)+abs(autocorr(goingdown(rightind)))/(abs(autocorr(goingdown(rightind)))+abs(autocorr(goingdown(rightind)+1)));
        width(jj)=rightpt-leftpt;
    %%%%%%
    %%%%%%
end
autocorrelations(kk,:)=median(autocorrstor);
widths(kk,:)=width;
[xpsd,psdense]=jcpsd2(pitches,8000);
psds(kk,:)=median(psdense);
end
g=8;
