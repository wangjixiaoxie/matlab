%This plotting program is used to zoom in on individual responses in a song
%file
%plot has three columns, leftmost column is just divided into a point, with
%labels for each of the four , center column is the song, right column is
%imagevec, with red lines, to indicate syllable boundaries.

inputbnds=[2.5 7.2]
numdiv=1;
song_fs=44100;
shftfs=1000;
binsize=.005
matfile='sum.mat'
load /cobain4/twarren4/g100o55/stim/stim.mat
load sum.mat
%stimnames={'p' 'm' 'm2' 'b' 'j' 'jc' 'j2' 'jc2'}
%stimnames={'p' 'm' 'b' 'j' 'jc'}
%load -mat '/cobain4/twarren/g100o55/stim/stim.mat' imagevec rawsong xlist randvec dat_shift
stimlist=[1 3  5 6 4];
stimnames={'m1.2' 'jig' 'jigc' 'p1.2'};
%stimnames={'jig' 'jigc' 'jig2' 'jig2c'}
clustnum=2;
rsptimes=epochsum(clustnum).rsplst;
meanvals=epochsum(clustnum).means(:,stimlist);
colors='bkkr'
linestyles={'-' '-' '--' '-'}

totleng=stimleng-2;
if(numdiv==2)
    %time_bnds=[2 (totleng+overlap)/numdiv+2 stimleng];
else
    time_bnds=inputbnds
end    
    rawsong=corpshiftednormg{5};

stimleng=length(rawsong)/song_fs;

plotrespfn2(rsptimes,meanvals,colors,rawsong,numdiv,stimleng,randvecfin,song_fs,stimnames,xlist,linestyles,time_bnds);
