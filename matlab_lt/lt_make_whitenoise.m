clear all; close all
duration =0.04; % in sec
filename='WN40_LT_Lo3000.wav'; % file saved in curr dir.

% PROBLEM - filtering - giving me bunch of NaNs/.

%% params
fs = 44050; % sound sample rate
% fs = 32000
    
%% MAKE
white=rand(ceil(fs*duration),1)*2-1;

%% Filter WN (optional);
lowFreq=100;
highFreq=3000;
PassOption = 3; % if 1, then bandstop filters WN, if 0, then bandpass filters, if 2, then WN is high, if 3, then Wn is low

Nq=ceil(fs/2); % nyquist freq

if PassOption==1; % bandstop filter
    % problem with this below: gives very high n. so that get a lot of
    % Nans in the filtered noise.  simpler to choose an n of 8.
%     Wp=[(lowFreq-200)/Nq (highFreq+200)/Nq]; % outer boundaries, 200 hz attenuation region.
%     Ws=[lowFreq/Nq highFreq/Nq]; % inner boundaries of stop band.
%     
%     [n,Wn]=buttord(Wp,Ws,3,60);
%     
%     [b,a]=butter(n,Wn,'stop');
%     
% %     freqz(b,a,512,fs);
% white_filt=filter(b,a,white);


  [b,a]=butter(8,[lowFreq*2/fs, highFreq*2/fs],'stop');
  white_filt=filtfilt(b, a, white);

end


% bandpass (i.e. WN only in a band)
if PassOption==0;
white_filt=bandpass(white,fs,lowFreq,highFreq,'hanningfir');
end


if PassOption==2; % high WN
white_filt=bandpass(white,fs,lowFreq,highFreq,'hanningfir');
end

if PassOption==3; % high WN
white_filt=bandpass(white,fs,lowFreq,highFreq,'hanningfir');
end

%% plot WN

% Freq domain data, unfiltered
y=white;

NFFT = 2^nextpow2(length(y));
spec=fft(y,NFFT)/length(y);
f=fs/2*linspace(0,1,NFFT/2+1);
PowSpec=spec.*conj(spec);

figure; plot(f,PowSpec(1:length(f)));
xlabel('Freq (hz)')
title('Power spectrum - unfiltered');


% filtered
try
   y=white_filt;

NFFT = 2^nextpow2(length(y));
spec=fft(y,NFFT)/length(y);
f=fs/2*linspace(0,1,NFFT/2+1);
PowSpec=spec.*conj(spec);

figure; plot(f,PowSpec(1:length(f)));
xlabel('Freq (hz)')
title('Power spectrum - filtered');
catch err % i.e. no filtered data exists.
end



% % time-domain data
% figure; plot(white)
% title('time domain');
% 
% % Plot filtered - optional
% spec=fft(white_filt);
% PowSpec=spec.*conj(spec);
% figure; plot(PowSpec)
% title('Filtered - power spectrum');

%% WRITE

% unfiltered WN
audiowrite(filename,white,fs);

% filtered WN
audiowrite(filename,white_filt,fs);

