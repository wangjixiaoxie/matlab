%lowpasswn.m
%commented 10/05/11 (SY)


%make wn file
%length of stimulus.
t=.06
totaltime=0.085
fs=44100
%lowpassfrequency
lopsvl=[10000];

stopvl=[4000 6500]
%tukey ratio
tukratio=.08;
type='low'
npts=fs*t
rawdata=rand(npts,1)*2-1
stop=0;
low=1;

%analyze wn file
rawdatashort=rawdata(1:1024);
fftout=fft(rawdatashort)
    fftfreqs=0:((44100/2)/(length(rawdatashort)/2)):44100/2;
    fftfreqs=fftfreqs(1:end-1)
    fftout=fftout(1:end/2)
    absout=abs(fftout.^2);
    figure
    plot(fftfreqs, (absout))
    
    
%now make notch file...THIS LINE 
[b,a]=cheby1(8,1,[lopsvl/(fs/2)],'low')
rawdataflt=filtfilt(b,a,rawdata);

if(stop==1)
    [b,a]=cheby1(8,1,[stopvl/(fs/2)],'stop')
    rawdataflt=filtfilt(b,a,rawdataflt);
end

    %normalize
    mxvl=max(abs(rawdataflt))
    rawdatafltnrm=rawdataflt/(mxvl*1.0001);
    lng=length(rawdataflt);
    tukout=tukeywin(lng,tukratio);
    rawdatafltnrm=rawdatafltnrm.*tukout;
    rawdatafltnrm=[rawdatafltnrm;zeros((totaltime-t)*fs,1)];
    rawdatafltnrm(1)=0;
    filename=['wn' num2str(t*1000) type num2str(lopsvl(1)) '.wav']
    
    wavwrite(rawdatafltnrm,44100,16,filename)
    figure; plot((1:totaltime*fs)/fs,rawdatafltnrm);
    
%  %now analyze wnfile played out over silence
%  cbinfile='r58pu68_111207_2022.2.cbin'
%             
%  [cbin,fs]=ReadCbinFile(cbinfile);
%  tstart=1.995
%  cbinshort=cbin(tstart*fs:tstart*fs+2047)
%  fftout=fft(cbinshort)
%     fftfreqs=0:((fs/2)/(length(cbinshort)/2)):fs/2;
%     fftfreqs=fftfreqs(1:end-1)
%     fftout=fftout(1:end/2)
%     absout=abs(fftout.^2);
%     figure
%     plot(fftfreqs, (absout))
%     
%     wncbinfile='r58pu68_111207_2037.167.cbin'
%  [wncbin,fs]=ReadCbinFile(wncbinfile);
%  tstart=5.96
%  wncbinshort=wncbin(tstart*fs:tstart*fs+511)
%  fftout=fft(wncbinshort)
%     fftfreqs=0:((fs/2)/(length(wncbinshort)/2)):fs/2;
%     fftfreqs=fftfreqs(1:end-1)
%     fftout=fftout(1:end/2)
%     absout=abs(fftout.^2);
%     figure
%     plot(fftfreqs, (absout))
%     
%     