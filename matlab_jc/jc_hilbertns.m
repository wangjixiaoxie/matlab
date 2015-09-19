function filtfreqvals=jc_hilbertns(shifted,lowfreq,highfreq,fs)
% This takes the hilbert transform
% all fluctuations faster than ms_period (in ms) are ignored

% It looks pretty crappy if you don't filter first.
sampling=fs; %44100;
Nyq=sampling/2;
lowratio=lowfreq/Nyq;
highratio=highfreq/Nyq;
%lowpassfilt=(1/ms_period)*1000*(1/Nyq);


[b,a]=butter(4,[lowratio highratio],'bandpass'); 
for i=1:size(shifted,1)
    x=shifted(i,:);
    y=filtfilt(b,a,x);
    % instantaneous deltapitch estimate
    idp=diff(unwrap(angle(hilbert(y))));
    % convert instantaneous deltapitch (radians/point) to frequency
    freqvals=sampling./((1./idp)*2*pi);
    %[c,d]=butter(4,lowpassfilt,'low');
    %filtfreqvals(:,i)=filter(c,d,freqvals);
    filtfreqvals(:,i)=freqvals;
end

    
    