% unzip chronux package into home/matlab/chronux
% add to path

% load ps9.mat

% Problem 1
% A. 
figure;plot(lfptimes(1:10000),lfpdata(1:10000))    % Figure1.png
figure;plot(lfptimes(1:1000),lfpdata(1:1000))    % Figure2.png
% Looks fairly oscillatory - period of ~ 100ms

% B.
% initialize 
    params={};
    params.Fs=1500; % sampling rate
    params.fpass=[0 150]; % low-pass filter - good for LFP analysis
    params.err=0;
    params.trialave=0;
% calculate spectrum - average power at different frequencies
    [S,f]=mtspectrumc(lfpdata(1:100000),params);
    % Note that I get error message if I use all the data because filters
    % of length greater than 2^20 (~10^6) are not supported.
    figure;plot(f,S)    % Figure3.png
    xlim([0 10])        % Figure4.png
    
% C.
    % Looks like interesting stuff at around 8Hz (period of ~120ms)
    
    
% Problem 2
% A.
% Spectrogram will show any changes in spectral structure over time.
    params.tapers=[]; % necessary to avoid error message
    [S t f]=mtspecgramc(lfpdata,[.5 .25],params);
    figure;imagesc(t,f,log(S'))     % Figure5.png

% B.
% rows are time points
% columns are frequencies
f(6:7)= [    7.3242    8.7891  ]; % Which freq bins to look at
thetapower=mean(S(:,6:7)'); % Theta power over time;
figure;plot(t,thetapower)       % Figure6.png
% It is necessary to calibrate times
    % speed is sampled at 30Hz
    % t is sampled at 4Hz
postimescal=postimes-min(postimes);
speedmagnitude=abs(speed);
% dumb interpolation
for i=1:length(t)
    tvals=find(postimescal>t(i)-0.125 & postimescal<t(i)+0.125);
    speedval(i)=mean(speedmagnitude(tvals));
end
figure;plot(speedval,thetapower,'*')    % Figure7.png
[r,p]=corrcoef(speedval,thetapower) = 0.0261, p=0.11
% Not a significant relationship

% Problem 3 - designing a filter with fdatool
% A. Astop1=60, Apass=1, Astop2=80
% These settings mean that parts of the signal between 6Hz and 12Hz are let
%   through 60 times as much as parts of the signal slower than 4Hz and 
%   80 times as much as parts of the signal faster than 14Hz
% ORDER = 1898

% B.
% The order decreases to 1063 if you don't need as much relative
% attenuation.  The calculation also becomes faster.

% C.
% Relative to the answer in B., the order doubles if you require the slope 
% to be twice as steep on both sides (i.e. attenuation of everything
% outside 5-13Hz instead of 4-14Hz as before).

% Problem 4
% A.
% Based on my version of matlab, the wording of the problem is backwards 
    % thetafilter.tf.num is B
    % thetafilter.tf.den is A
    % Because filter(B,A,X) is the way the matlab command works
    
lfpFILT=filter(thetafilter.tf.num,thetafilter.tf.den,lfpdata);
figure;hold on; % figure8.png
plot(lfptimes,lfpdata)
plot(lfptimes,lfpFILT,'r','Linewidth',2)
xlim([min(lfptimes) min(lfptimes)+2]) % look at first 2 seconds
% looks bad

figure;hold on; % figure9.png
plot(lfptimes,lfpdata)
plot(lfptimes,lfpFILT,'r','Linewidth',2)
xlim([min(lfptimes)+30 min(lfptimes)+32]) % look at 30-32 seconds
% looks better
% Edges reflect properties of the filter instead of properties of the
% signal.  The filter is 1273 points wide - the left half of the filter is 
% about 0.42 seconds wide.

% B.
lfpFILTFILT=filtfilt(thetafilter.tf.num,thetafilter.tf.den,lfpdata);
figure;hold on; % figure10.png
plot(lfptimes,lfpdata)
plot(lfptimes,lfpFILTFILT,'r','Linewidth',2)
xlim([min(lfptimes) min(lfptimes)+2]) % look at first 2 seconds
% This does a much better job at the edges.

% Problem 5

hilbertdata=hilbert(lfpdata);
magnitude_hilbert=abs(hilbertdata);
phase_hilbert=angle(hilbertdata);

% A
figure;hold on;
subplot(211);hold on; % figure11.png
plot(lfptimes,abs(lfpdata))
plot(lfptimes,magnitude_hilbert,'r')
xlim([min(lfptimes)+4 min(lfptimes)+5])
% Note that magnitude estimated by hilbert is often greater than actual
% magnitude of the signal because some of the energy is transferred into 
% phase energy.

% B
subplot(212);hold on
plot(lfptimes,phase_hilbert,'g')
plot(lfptimes,lfpdata/300)
plot([min(lfptimes)+4 min(lfptimes+5)],[0 0],'k')
xlim([min(lfptimes)+4 min(lfptimes)+5])

% C
% Look up times
for i=1:length(spiketimes)
    [b,lfpind(i)]=min(abs(lfptimes-spiketimes(i)));
end
figure;hist(phase_hilbert(lfpind)) % Figure12.png
% Spikes tend to happen near -pi and pi 
% There is a tendency for spikes occur at peaks of theta rhythm.

% D
fs=1500;
raweeg=zeros(length(spiketimes),fs*2+1);
for i=1:length(spiketimes)
    raweeg(i,:)=lfpdata(lfpind(i)-fs:lfpind(i)+fs);
end
relativetime=[-fs:fs]/1500;
figure;plot(relativetime,mean(raweeg)) % figures 13 and 14
% Spikes tend to happen at the negative peaks of the LFP

% E
for i=1:length(spiketimes)
    [b,posind(i)]=min(abs(postimes-spiketimes(i)));
end
positions=pos(posind);
thetaphases=phase_hilbert(lfpind); % from C.
figure;plot(positions,thetaphases,'*') % figure 15

% F
count=0;
% Positions with + velocity
diffs=diff(pos);
% Positions with a spike that have + velocity
increasingindex=diffs(posind)>0;
figure;plot(posind(increasingindex),thetaphases(increasingindex),'*') % Fig 16
figure;plot(posind(increasingindex),thetaphases(increasingindex)+2*pi,'*') % Fig 17



