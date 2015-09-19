%make wn file
t=.08
fs=44100
npts=fs*t
rawdata=rand(npts,1)*2-1
wavwrite(rawdata,44100,16,'wn60.wav')

%analyze wn file
rawdatashort=rawdata(1:2048);
fftout=fft(rawdatashort)
    fftfreqs=0:((44100/2)/(length(rawdatashort)/2)):44100/2;
    fftfreqs=fftfreqs(1:end-1)
    fftout=fftout(1:end/2)
    absout=abs(fftout.^2);
    figure
    plot(fftfreqs, (absout))
    
    
%now make notch file
[b,a]=cheby1(8,1,[6500*2/fs 8500*2/fs],'stop')
rawdataflt=filtfilt(b,a,rawdata);
rawdatashortflt=rawdataflt(1:2048)
fftout=fft(rawdatashortflt)
    fftfreqs=0:((44100/2)/(length(rawdatashortflt)/2)):44100/2;
    fftfreqs=fftfreqs(1:end-1)
    fftout=fftout(1:end/2)
    absout=abs(fftout.^2);
    figure
    plot(fftfreqs, (absout))
    
    %now write this wav file
    
    %normalize
    mxvl=max(abs(rawdataflt))
    rawdatafltnrm=rawdataflt/mxvl;
    wavwrite(rawdatafltnrm,44100,16,'wn80filt.wav')
    
    
 %now analyze wnfile played out over silence
 cbinfile='r58pu68_111207_2022.2.cbin'
            
 [cbin,fs]=ReadCbinFile(cbinfile);
 tstart=1.995
 cbinshort=cbin(tstart*fs:tstart*fs+2047)
 fftout=fft(cbinshort)
    fftfreqs=0:((fs/2)/(length(cbinshort)/2)):fs/2;
    fftfreqs=fftfreqs(1:end-1)
    fftout=fftout(1:end/2)
    absout=abs(fftout.^2);
    figure
    plot(fftfreqs, (absout))
    
    wncbinfile='r58pu68_111207_2037.167.cbin'
 [wncbin,fs]=ReadCbinFile(wncbinfile);
 tstart=5.96
 wncbinshort=wncbin(tstart*fs:tstart*fs+511)
 fftout=fft(wncbinshort)
    fftfreqs=0:((fs/2)/(length(wncbinshort)/2)):fs/2;
    fftfreqs=fftfreqs(1:end-1)
    fftout=fftout(1:end/2)
    absout=abs(fftout.^2);
    figure
    plot(fftfreqs, (absout))
    
    