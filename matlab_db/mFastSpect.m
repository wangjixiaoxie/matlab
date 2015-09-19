function [spec,t,f]=mFastSpect(rawsong,Fs,SPTH,nfft,olap,sm_win,F_low,F_High);
% [smooth,spec,t,f]=evsmooth(rawsong,Fs,SPTH,nfft,olap,sm_win,F_low,F_High);
% returns the smoothed waveform/envelope + the spectrum
%

if (~exist('F_low'))
    F_low  = 750.0;
end
if (~exist('F_high'))
    F_high = 15000.0;
end
if (~exist('nfft'))
    nfft = 512;
end
if (~exist('olap'))
    olap = 0.8;
end
if (~exist('sm_win'))
    sm_win = 2.0;%ms
end
if(~exist('SPTH'))
    SPTH=0.01;
end

filter_type = 'hanningffir';

filtsong=bandpass(rawsong,Fs,F_low,F_high,filter_type);

if (nargout > 1)
    spect_win = nfft;
    noverlap = floor(olap*nfft);
    [spec,f,t] = specgram(filtsong,nfft,Fs,spect_win,noverlap);
    if (exist('SPTH'))
        if(SPTH == 0)
            SPTH = mean(mean(abs(spec)));
            p = find(abs(spec)<=SPTH);
            spec(p) = SPTH / 2;
        else
            p = find(abs(spec)<=SPTH);
            spec(p) = SPTH;
        end
    end
end


return;
