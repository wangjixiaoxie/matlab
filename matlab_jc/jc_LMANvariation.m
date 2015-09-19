function LMANcurves=jc_LMANvariation(x1,x2)
% e.g. x1=ACSF(3).residHILB2(1900:2500,21:50);
% e.g. x2=INA(3).residHILB2(2000:2600,1:30);

%%%%%%%%
% The idea of this is two-fold.
% 1. produce the autocorrelation of the "LMAN-dependent variability"
% 2. simulate LMAN variability -- decompose signals containing and lacking
% LMAN, subtract the spectral difference and add that back to the signal
% containing LMAN --- TO BE USED for learning simulation




%%%%% 12-15-08 %%%%%%%%
    % ********normalize after taking the mean so that you take into account the
    % importance of the notes with big residuals
figure;plot(mean(jc_xcorr(x1)')./max(mean(jc_xcorr(x1)')),'r')
hold on;plot(mean(jc_xcorr(x2)')./max(mean(jc_xcorr(x2)')),'b')
plotLMANxcorr(x1,x2);

% This code plots the widest it can be
for i=1:size(x1,1)
    ll(i)=10;
end
hold on;plot(xcorr(ll)./max(xcorr(ll)),'k')

%%% Basic method 
    % This takes noLMAN (x2) and adds randomized magnitude differences between 
    % yesLMAN (x1) and noLMAN (x1-x2) at randomized phases selected from the
    % noLMAN phase distribution (why this distn?).  This sum should approximately mimic the
    % timescale of yesLMAN variation, as confirmed by similarity in the
    % autocorrelation function.
long=min([size(x1,2) size(x2,2)]);
for i=1:long
    a=x1(:,i); 
    b=x2(:,i); 
    z1=fft(a);
    z2=fft(b);
    mag1(i,:)=abs(z1);
    mag2(i,:)=abs(z2);
    pha1(i,:)=angle(z1);
    pha2(i,:)=angle(z2);
end
mag=mag1-mag2; %%%%%%
%Simulation;
count=0;
for j=1:10
for i=1:100
    count=count+1;
    random=round(rand*(long-1))+1;
    m=2*mag(round(rand*(long-1))+1,:);
    p=pha2(round(rand*(long-1))+1,:); % this could be changed.
    re=m.*cos(p);     
    im=m.*sin(p);
    dat=complex(re,im);
    %figure;plot(real(ifft(dat)))
    gg=real(ifft(dat));
    gg1(:,i)=x2(:,random)+gg';
    acorrest(:,i)=xcorr(gg1(:,i));
    LMANcurves(count,:)=gg;
end
figest(j,:)=mean(acorrest')./max(mean(acorrest'));
end

hold on;plot(figest','g')
