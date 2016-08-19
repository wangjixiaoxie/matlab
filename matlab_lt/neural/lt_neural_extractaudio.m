%% 7/21/16 - LT extract analog audio from Intan recording

function lt_neural_extractaudio(filename)

%% ==========

% filename='MultBirds_AudioOnly_160721_175929.rhd';

% put all files in folder into batch
!ls *.rhd* > batch
batch='batch';


%% === single file
figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

fline=filename;
    
    % ==========
    [amplifier_data,board_dig_in_data,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, board_dig_in_channels] = pj_readIntanNoGui(fline);
    
    fs=frequency_parameters.board_adc_sample_rate;
    % --- NOTE ON THE OUTPUT VARIABLES
    
    
    % ---
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(fline);
    tt=[1:length(board_adc_data)]/fs;
    plot(tt, board_adc_data, 'k');
    axis('tight');
    
    % === get params
    window = 0.016*fs; % make it 16ms, to match evtaf stuff
    nfft=2^nextpow2(window); % go to next pow 2 size
    olap=0.8;
    noverlap=olap*window;
    
    
    % ==== filter song
    songdat=board_adc_data;
	filter_type='hanningfir';
	F_low  = 500;
	F_high = 8000;
    songdat=bandpass(songdat,fs,F_low,F_high,filter_type);
    
    
    % === colect spectrogram
    [sp, f, t] = spectrogram(songdat, window, noverlap, nfft, fs);
    sp=abs(sp);
    
    % ========= cut off lower higher frequencies
    maxfreq=7500;
    inds=f<maxfreq;
    f=f(inds);
    sp=sp(inds, :);
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(fline);
       
    
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
    if (0)
                   sptemp_vector=reshape(sptemp, numel(sptemp), 1);
            Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
            [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
            plot(Xcenters, Ybinned);
            
            [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
            Cmin=Xcenters(maxind);
            Cmax=max(sptemp_vector);
            
            % Plot
            colormap('hot')
            imagesc(t, f, sptemp, [Cmin+Cmin/10 Cmax]);        
    end
    
    % Plot
    imagesc(t, f, sptemp);
    
%     axis([t(1) t(end) f(end) f(1)]);
    axis([t(1) t(end) f(1) f(end)]);
    
%     set(gca, 'YTick', linspace(f(1), f(end), 5));
%     set(gca, 'YTickLabel', linspace(f(end), f(1), 5));
    

%% 
fid=fopen(batch);

fline=fgetl(fid);

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


while ischar(fline)
    
    % ==========
    [amplifier_data,board_dig_in_data,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, board_dig_in_channels] = pj_readIntanNoGui(fline);
    
    % --- NOTE ON THE OUTPUT VARIABLES
    
    
    % ---
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(fline);
    plot(board_adc_data, 'k');
    axis('tight');
    
    % === get params
    fs=frequency_parameters.board_adc_sample_rate;
    window = 0.016*fs; % make it 16ms, to match evtaf stuff
    nfft=2^nextpow2(window); % go to next pow 2 size
    olap=0.8;
    noverlap=olap*window;
    
    
    % ==== filter song
    songdat=board_adc_data;
	filter_type='hanningfir';
	F_low  = 500;
	F_high = 8000;
    songdat=bandpass(songdat,fs,F_low,F_high,filter_type);
    
    
    % === colect spectrogram
    [sp, f, t] = spectrogram(songdat, window, noverlap, nfft, fs);
    sp=abs(sp);
    
    % ========= cut off lower higher frequencies
    maxfreq=7500;
    inds=f<maxfreq;
    f=f(inds);
    sp=sp(inds, :);
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(fline);
       
    
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
    if (0)
                   sptemp_vector=reshape(sptemp, numel(sptemp), 1);
            Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
            [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
            plot(Xcenters, Ybinned);
            
            [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
            Cmin=Xcenters(maxind);
            Cmax=max(sptemp_vector);
            
            % Plot
            colormap('hot')
            imagesc(t, f, sptemp, [Cmin+Cmin/10 Cmax]);        
    end
    
    % Plot
    imagesc(t, f, sptemp);
    
%     axis([t(1) t(end) f(end) f(1)]);
    axis([t(1) t(end) f(1) f(end)]);
    
%     set(gca, 'YTick', linspace(f(1), f(end), 5));
%     set(gca, 'YTickLabel', linspace(f(end), f(1), 5));
    
    % ============ get next file
    fline=fgetl(fid);
end
