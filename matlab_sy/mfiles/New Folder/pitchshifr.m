%First part prepares the aif file.

infile='r82pk21_031005_2049.38.cbin'
mintime=7;
maxtime=10;
%pitchshifts
sr=44150.1103752759372000

FFS=44100
[dat]=evsoundin('',infile,'obs0');

%adjust the data, so that it's sampled at 44.1 khz
redat=interp1([0:length(dat)-1]/sr,dat,[0:(1/FFS):((length(dat)-1)/sr)],'linear');
tt=[1:length(redat)]/FFS;
pp=find(tt>mintime&tt<maxtime);

redat=redat(pp);

displaybuffer=500;

%Here I normalize song to write aif file
datnorm=redat/(1.01*max(abs(dat)));
%This is requirement of aif converter
dataif=32768*datnorm;
aiffwrite('song.aif',dataif,FFS,16);


%The following code will read the aif file back in.
%pitchshifts=[6 12 .6 .12]
%for i=1:length(pitchshifts)
    datsame=aiffread('song.aif')';
    [datp2,fs]=aiffread('songp2.aif');
    datp10=aiffread('songp10.aif');
    datp50=aiffread('songp50.aif');
    datm2=aiffread('songm2.aif');
    datm10=aiffread('songm10.aif');
    datm50=aiffread('songm50.aif');
    datp50m33=aiffread('songp50m33.aif');
    %/*now I can pitchshift*/



%By adding more points with the sampling rate
syllint=dat
syllint=interp(dat(ppp),intlevel)

%now it is too long, so I am going to cut off and bottom parts.









%By taking away points, with the same sampling rate
