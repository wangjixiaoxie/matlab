function lt_plot_spectrogram(songdat, fs, bk_red, plot_ms)
%% lt 8/18/16 - takes song vector and plots spectrogram on current active axes

% plot_ms=1, ms (if 0, sec)
    
if isempty('plot_ms')
    plot_ms=0;
end

%% params
% songdat; % vector
% fs, samp,sec
% bk_red = plots black and red

if isempty('bk_red')
    bk_red=0;
end
   

% - spectrogram
window = 0.016*fs; % make it 16ms, to match evtaf stuff
nfft=2^nextpow2(window); % go to next pow 2 size
olap=0.8;
noverlap=ceil(olap*window);

% - filter
filter_type='hanningfirff';
F_low  = 500;
F_high = 8000;

%% run

% ==== filter song
songdat=bandpass(songdat,fs,F_low,F_high,filter_type);

% === colect spectrogram
[sp, f, t] = spectrogram(songdat', window, noverlap, nfft, fs);
if plot_ms==1
t=t.*1000;
end
sp=abs(sp);

% ========= cut off lower higher frequencies
maxfreq=7500;
inds=f<maxfreq;
f=f(inds);
sp=sp(inds, :);

% TAKE LOG of sp
% first, convert any sp values of 0 to non-zero(to the lowest value present);
% solves problem of taking log of 0
pp=find(sp>0);
mntmp = min(min(sp(pp)));
pp=find(sp==0);
sp(pp) = mntmp;

% second, take log
sptemp=log(sp);
sptemp = sptemp - min(min(sptemp));
sptemp = ((2^8) - 1)*(sptemp./max(max(sptemp)));

% === black and red
if bk_red==1
    sptemp_vector=reshape(sptemp, numel(sptemp), 1);
    Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
    [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
%     plot(Xcenters, Ybinned);
    
    [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
    Cmin=Xcenters(maxind);
    Cmax=max(sptemp_vector);
    
    % Plot
    colormap('hot')
    imagesc(t, f, flipud(sptemp), [Cmin+Cmin/10 Cmax]);
else
            imagesc(t, f, sptemp);
end
            axis([t(1) t(end) f(1) f(end)]);
