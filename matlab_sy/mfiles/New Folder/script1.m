intlevel=2


[dat]=evsoundin('','r82pk21_031005_2049.38.cbin','obs0');
pp=findstr(labels,'a');

displaybuffer=500;

tt=[1:length(dat)]/fs;

%for i=1:length(pp)

ppp=find((tt>=onsets(13)*1e-3)&(tt<=offsets(13)*1e-3));


%what is this 100 referring to.
%[sm,sp,t,f]=evsmooth(dat(ppp(1)-500:max(ppp)+500),fs,100);
figure
[sm,sp,t,f]=evsmooth(dat_6,fs,100);

imagesc(t,f,log(abs(sp)));
set(gca,'YD','n')
colorbar
%
wavwrite(dat,fs,'song.wav')
datrsmple=resample(dat,44150000,44150110);
dat=32768*dat;

aiffwrite('song.aif',dat,fs,16);
dat_6=aiffread('song-out.aif');
resample
aiffread

/*now I can pitchshift*/



%By adding more points with the sampling rate
syllint=dat
syllint=interp(dat(ppp),intlevel)

%now it is too long, so I am going to cut off and bottom parts.









%By taking away points, with the same sampling rate
