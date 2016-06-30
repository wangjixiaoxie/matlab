%% ===== taking signal, decomposing by fft, etc.

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% ==== make signal

N=2^8;


x1=(1-2*rand(1, N))/2;
k=5; x2=2*sin(2*pi*k*(1:N)/N);
k=2; x3=5*sin(2*pi*k*(1:N)/N);
k=71; x4=4*sin(2*pi*k*(1:N)/N);

X=x1+x2+x3+x4;


% ==== plot signal
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('raw signal');
plot(1:N, X, '.-')


%% hamming window?
% if (0)
ham=hamming(N);
X=X.*ham';
% end

%% perform FFT

NFFT=2^10;
% NFFT=2^8;

Y=fft(X, NFFT);



%% ==== plot sin and cos components

x=linspace(0,1,NFFT/2+1)*N/2; % number of samples is based on NFFT. ff bins are based on number of samples in original data (here N is effectively Fs, giving
% us frequencies that are in units of cylces/sec.

YY=Y(1:NFFT/2+1);


Ycos=real(YY);
Ysin=imag(YY);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('cosine component')
plot(x,Ycos, '.-');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('sine component')
plot(x,Ysin, '.-');

% === plot single sided amplitude
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('single sided amplitude spectrum');

YY=abs(Y(1:NFFT/2+1));
plot(x,YY, '.-k')

% === plot single sided power
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('single sided power');

YY=abs(Y(1:NFFT/2+1));
plot(x,YY.^2, '.-k')



%% === example on doc
% 
% Fs = 1000;                    % Sampling frequency
% T = 1/Fs;                     % Sample time
% L = 1000;                     % Length of signal
% t = (0:L-1)*T;                % Time vector
% % Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
% x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
% y = x + 2*randn(size(t));     % Sinusoids plus noise
% plot(Fs*t(1:50),y(1:50))
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('time (milliseconds)')
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
