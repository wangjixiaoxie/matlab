
%make wn file
%length of stimulus.
t=.05
fs=44100
npts=fs*t
rawdataflt=rand(npts,1)*2-1


    %normalize
    mxvl=max(abs(rawdataflt))
    rawdatafltnrm=rawdataflt/(mxvl*1.0001);
    lng=length(rawdataflt);
    rawdatafltnrm(1)=0;
    filename=['wn' num2str(t*1000)  '.wav']
    
    wavwrite(rawdatafltnrm,44100,16,filename)
    
    