%lowpasswn.m
%commented 3.27.08


%make wn file
%length of stimulus.
t=.05
fs=44100
%lowpassfrequency
lopsvl=4500;
%tukey ratio
tukratio=.02;

npts=fs*t
rawdata=rand(npts,1)*2-1


% %analyze wn file
% rawdatashort=rawdata(1:2048);
% fftout=fft(rawdatashort)
%     fftfreqs=0:((44100/2)/(length(rawdatashort)/2)):44100/2;
%     fftfreqs=fftfreqs(1:end-1)
%     fftout=fftout(1:end/2)
%     absout=abs(fftout.^2);
%     figure
%     plot(fftfreqs, (absout))
    
    
%now make notch file
[b,a]=cheby1(8,1,[lopsvl],'low')
rawdataflt=filtfilt(b,a,rawdata);

    %normalize
    mxvl=max(abs(rawdataflt))
    rawdatafltnrm=rawdataflt/(mxvl*1.0001);
    lng=length(rawdataflt);
    tukout=tukeywin(lng,tukratio);
    filename=['wn' num2str(t*1000) 'sub' num2str(lopsvl) '.wav']
    
    wavwrite(rawdatafltnrm,fs,16,'wn80filt.wav')
    
    
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