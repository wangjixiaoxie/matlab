function [smooth] = mquicksmooth(rawsong,Fs,DOFILT,sm_win,F_low,F_high)
%
% [smooth]=mquicksmooth(rawsong,Fs,DOFILT,,sm_win,F_low,F_High);
% returns the smoothed waveform/envelope 
%
% Like evsmooth() but doesn't compute or return spectrogram, for speed
% purposes.
%
%
%


filter_type='hanningfir';

if (~exist('F_low','var'))
	F_low  = 500.0;
end
if (~exist('F_high','var'))
	F_high = 10000.0;
end
if (~exist('nfft','var'))
	nfft = 512;
end
if (~exist('sm_win','var'))
	sm_win = 2.0;%ms
end

if (~exist('DOFILT','var'))
    DOFILT=1;
end

if (DOFILT==1)
    filtsong=bandpass(rawsong,Fs,F_low,F_high,filter_type);
else
    filtsong=rawsong;
end
    

squared_song = filtsong.^2;

%smooth the rectified song
len = round(Fs*sm_win/1000);
h   = ones(1,len)/len;
smooth = conv(h, squared_song);
offset = round((length(smooth)-length(filtsong))/2);
smooth=smooth(1+offset:length(filtsong)+offset);
return;
