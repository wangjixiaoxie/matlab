function [x]=jc_fluc(ffreq1)
ffreq1=200;
fs = 32000; % cbin sampling rate
t = 0:1/fs:0.03125; %3200pts - 100ms
modu=0.2*cos(ffreq1*pi*t);
% est=modu*60+2500;
y=0.2*vco(pp,[2350 2400],fs);
x = 0.2*vco(modu,[2350 2400],fs); % 2400 and 3600 are the minimum and maximum frequency values that the pitch oscillates between
xfinal=[zeros(1,512) x zeros(1,512)];
pitch=jc_pitchmat1024(xfinal,1024,1020,1,2200,2900,[1],'obs0',1);

xfinal=[zeros(1,512) x zeros(1,512)];
rawDIRfinal=[zeros(1,512) rawchunk(2,1:1000) zeros(1,512)];
    



%z= 0.8*vco(cos(ffreq2*pi*t),[2940 3060],fs);
%k=3000+60*(0.2*cos(ffreq1*pi*t)+0.8*cos(ffreq2*pi*t));
%y=noise*rand(size(x));  % add in some noise at each point in time
x=x+z;
figure; spectrogram(x,128,124,128,32000, 'yaxis')

%Gaussian windows (sonogram)
%sigma=2;
%OVERLAP=1010;
%N=1024;
%t=-N/2+1:N/2;
%sigma=(sigma/1000)*fs;
%w=exp(-(t/sigma).^2);
%q=specgram(x,1024,[],w,OVERLAP)+eps;
%sonogram=abs(flipdim(q,1));