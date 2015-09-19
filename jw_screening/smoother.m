function smoothed=smoother(syl,Fs,sm_win)
%smoothed=smoother(syl,Fs,sm_win)
%sm_win is in ms

filtsyl=bandpass(syl,Fs,300,8000);

%calculate square of signal: proportional to power
squared_song = filtsyl.^2;
  
%smooth the rectified song
len=round(Fs*sm_win/1000);                      
h=hann(len)/(.5*len);
smooth=conv(h, squared_song);
offset=round((length(smooth)-length(filtsyl))/2); %get rid of convolution induced offset
smoothed=smooth(1+offset:length(filtsyl)+offset);
