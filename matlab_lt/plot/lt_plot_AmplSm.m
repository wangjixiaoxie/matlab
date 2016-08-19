%% lt 8/17/16 - plot smooth ampl waveform, a la evsonganaly

function [smooth, tt]=lt_plot_AmplSm(dat, fs)

% dat is raw (pre filt)
% tt in sec

% NOTE: does not plot, just gives output
%%
	filter_type='hanningfir';
	F_low  = 500.0;
	F_high = 10000.0;
	sm_win = 2.0; %ms

filtsong=bandpass(dat,fs,F_low,F_high,filter_type);
squared_song = filtsong.^2;

%smooth the rectified song
len = round(fs*sm_win/1000);
h   = ones(1,len)/len;
smooth = conv(h, squared_song);
offset = round((length(smooth)-length(filtsong))/2);
smooth=smooth(1+offset:length(filtsong)+offset);

tt=(1:length(smooth))/fs;


