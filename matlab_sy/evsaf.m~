
%This is a rough script created approx. 1/2007, commented 1/2009
%by Tim Warren.  Purpose is to simulate parameters collected during evsaf2.

%Frequency bins.
Fs=32000
F_low  = 750.0;
F_high = 10000.0;
nfft = 256;

%Where are these used?
olap = 0.8;
noverlap = floor(olap*nfft);
filter_type = 'hanningfir';
startind=1
indlst=[];
keepdata=[];
ratlist=[];
entlist=[];
songindlst=[];
triglist=[];
onsets=0;
%test for trigger
trigtest=1

%Volume bnds
minvol=1000
maxvol=18000;

%entropy bounds.
maxent=5.5;
maxsamrat=.5
counter=3;
USEX=0;
CS='obs0'
%make a bunch of pts to start data.

%next step is to read in onsets and offsets
%and classify every bin as being in song  or out of song....
%from the cbin.not.mat, load the file onsets and offsets

%three batch files of cbins.

noisebatchfile='batchnoise';
songbatchfile='batchsong';
callbatchfile='batchcall';

bt=songbatchfile;

ff=load_batchf(bt);
for ii=1:length(ff)
	trigcount=0;
    fn=ff(ii).name
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1));
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf(fn,USEX);
	load([fn,'.not.mat']);

	[dat,fs]=evsoundin('',fn,CS);

    datfile=dat;

    startind=1
    indlst=[];
    keepdata=[];
    ratlist=[];
    entlist=[];
    vol_list=[];
    %this loop is for each datfile
    while(startind < length(datfile)-nfft)
        indlst=[indlst startind];
        startind=startind+(nfft)
    end 
    for jj=1:length(indlst)
        indlst(jj)
        datvals=datfile(indlst(jj):indlst(jj)+nfft-1);
        startind=startind+nfft;
        %filtsong=bandpass(datvals,Fs,F_low,F_high,filter_type);
        spect_win = nfft;
        noverlap = floor(olap*nfft);
%         [spec,f,t,p] = spectrogram(datvals,spect_win,noverlap,nfft,Fs);
    
        
        
		fdat=(abs(fft(hamming(nfft).*datvals))./nfft);
        ffv=get_fft_freqs(nfft,fs);
		ffv=ffv(1:end/2);
		frqind=find((ffv>=F_low)&(ffv<=F_high));

        
        
        %volume calculation
       
        
        volcalc=sum(fdat(frqind));
        vol_list=[vol_list volcalc]; 
        
        %entropy calculation
        pnorm=(fdat(frqind))/sum(fdat(frqind));
        entcalc=-(pnorm'*log2(pnorm));
        entlist=[entlist entcalc];
        
        %ratio calculation
        hi_F_id=find(ffv>10000);
        lo_F_id=find(ffv<=10000);
        mean_p_hi=mean(mean(fdat(hi_F_id,:))); 
        mean_p_lo=mean(mean(fdat(lo_F_id,:)));
        hi_over_lo=mean_p_hi/mean_p_lo;
   
        ratlist=[ratlist hi_over_lo];
       if (trigtest)
            
            if (hi_over_lo<maxsamrat)&(entcalc<maxent)&(volcalc<maxvol)&(volcalc>minvol)
                trigcount=trigcount+1;
            else
                trigcount=0;
            end
            if trigcount>(counter-1)
                triglist=[triglist jj];
            end
        end
 
        
    end
    songindlst=zeros(length(indlst),1);    
    if (onsets)
        
        for kk=1:length(onsets)
                ind=find((indlst/fs*1000)>(onsets(kk)+16)&(indlst/fs*1000)<(offsets(kk)-16));
                songindlst(ind)=1;
            end    
        end
    
    cmd = ['save ',fname,'.not.mat entlist ratlist indlst vol_list songindlst triglist -append']
        eval(cmd)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of main program



%plotting
fname='o26bk80_170407_0702.1706.cbin'
cmd=['load ' fname '.not.mat'];
eval(cmd);
[datfile,fs]=evsoundin('',fname,'obs0');

songind=find(songindlst==1);
nosongind=find(songindlst==0);


        

figure
ax2(1)=subplot(5,1,1)
 [sm,sp,t,f]=evsmooth(datfile,fs,.005);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);

ax2(2)=subplot(5,1,2)
% ind=find(ratlist>.1)
% rastlist=indlst(ind);
% rastlist=rastlist'
% rastlist(:,2)=1;
% plotrasters(rastlist);
    plot(indlst(songind)/fs,ratlist(songind),'b.');
    hold on;
    plot(indlst(nosongind)/fs,ratlist(nosongind),'r.');
    title('samratio')
ax2(3)=subplot(5,1,3);
    plot(indlst(songind)/fs,vol_list(songind),'b.');
    hold on;
     plot(indlst(nosongind)/fs,vol_list(nosongind),'r.');
    title('volume')

ax2(4)=subplot(5,1,4);

plot(indlst(songind)/fs,entlist(songind),'b.')
hold on;
plot(indlst(nosongind)/fs,entlist(nosongind),'r.')

    title('entropy')
ax2(5)=subplot(5,1,5);
    triglist2=indlst(triglist)'
   
    triglist2(:,2)=1;
    plotrasters(triglist2);
linkaxes(ax2,'x');
%This bit of secondary analysis will determine if a given calculated value
%is within song or not.  A (stupid) way to do this is to loop through
%onsets for each file...and then mark those points with a one if they fall
%in between onsets and offsets.

entlstsong=[];
ratlstsong=[];
volstsong=[];
entlstnoise=[];
ratlstnoise=[];
volstnoise=[];
for ii=1:length(ff)
	fn=ff(ii).name
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1));
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf(fn,USEX);
	load([fn,'.not.mat']);
    
    
    songindlst=zeros(indlst,1);
    for ii=1:length(onsets)
        ind=find((indlst/fs*1000)>(onsets(ii)+16)&(indlst/fs*1000)<(offsets(ii)-16))
        songindlst(ind)=1;
        
        
    end
    %now combine the song and noise values of the entlist and the ratlist
    songind=find(songindlst==1)
    noiseind=find(songindlst==0);
    
    entlstsong=[entlstsong entlist(songind)];
    ratlstsong=[ratlstsong ratlist(songind)];
    volstsong=[volstsong vol_list(songind)];
    entlstnoise=[entlstnoise entlist(noiseind)];    
    ratlstnoise=[ratlstnoise ratlist(noiseind)];
    volstnoise=[volstnoise vol_list(noiseind)];

end


%divide song into 512 ms chunks
    spect_win=32;
    nfft=round(fs*spect_win/1000)
    Fs=32000
    noverlap=0

    [spec,f,t] = spectrogram(datvals,nfft,Fs,spect_win,noverlap);
    p = find(abs(spec)<=SPTH);
    spec(p) = SPTH;
    end
    
    
    
    
    %where does this trigger during song
    
%Amplitude measurement.

%calculation of spectral entropy.