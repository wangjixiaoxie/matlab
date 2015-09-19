function [datfile]= mkwnfile(time,fs,freqbnd,outfilename)

    datfile=rand(time*fs,1);
    datfile=2*datfile;
    datfile=datfile-1;
    if freqbnd
       % [b,a]=cheby1(5,0.5,[freqbnd/(fs/2) ],'low')
        %datfile=filtfilt(b,a,datfile);
        %[b,a]=butter(2,freqbnd(2)/(fs/2),'low')
        %datfile=filtfilt(b,a,datfile);
    end
    datfile=.999*(datfile/max(abs(datfile)));
    datfile(1)=0
    datfile(end)=0
    wavwrite(datfile,fs,16,outfilename);