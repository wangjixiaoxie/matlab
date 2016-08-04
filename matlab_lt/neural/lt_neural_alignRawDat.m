%% LT 8/3/16 - single file, aligns raw signal for all channels desired by user
function lt_neural_alignRawDat

clear all;

DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
AnalogChans_zero=[0]; % assumes that this is audio
AmpChans_zero=[9 14 19 23];

neuralFiltLow=300;
neuralFiltHi=[];


filename='bk7_chan9_14_19_23_Audio_160803_170525.rhd';

%% === extract for this file

[amplifier_data,board_dig_in_data,frequency_parameters, board_adc_data, board_adc_channels, ...
    amplifier_channels, board_dig_in_channels] = pj_readIntanNoGui(filename);


fs=frequency_parameters.amplifier_sample_rate; % for neural, analog, and dig


%% ==== make neural filter
% butterworth bandpath
N=4;

if neuralFiltLow<min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency])
    % make low freq at least as high at acquisition frequency.
    neuralFiltLow=min([frequency_parameters.actual_lower_bandwidth, ...
        frequency_parameters.actual_dsp_cutoff_frequency]);
end

neuralFiltHi=frequency_parameters.actual_upper_bandwidth;

[filt_b,filt_a]=butter(N,[neuralFiltLow*2/fs, neuralFiltHi*2/fs]);


%% ---- get timebins
if exist('amplifier_data', 'var')
    tt=1:length(amplifier_data(1, :));
    tt=tt./fs; % in sec
else
    tt=1:length(board_adc_data(1, :));
    tt=tt./fs; % in sec
end
tt=tt*1000; % to ms;


%% ==== For each desired signal type, collect all waveforms
% will override input channel IDs and extracted data - so only keeping
% desired channels, and in correct order as in the input. If no input, then
% will keep all channels that exist.

% will also PLOT (raw, not filtered)


% ==== DIG CHANNELS




% ==== ANALOG (board aux) inputs
if isempty(AnalogChans_zero)
    
else
    chanstmp=[];
    dattmp=[];
    
    for i=AnalogChans_zero
        
        ind=find([board_adc_channels.chip_channel]==i); % ind to find this chan
        
        if isempty(ind)
            % then did not recoprd this chan...
            disp(['PROBLEM - do not have data for analog chan (skipping) ' num2str(i)]);
            continue
        end
        
        chanstmp=[chanstmp i];
        dattmp=[dattmp; board_adc_data(ind, :)];
    end
    
    AnalogChans_zero=chanstmp;
    board_adc_data=dattmp;
end



% ==== AMPLIFIER CHANNELS (neural)
if isempty(AmpChans_zero)
    
else
    chanstmp=[];
    dattmp=[];
    
    for i=AmpChans_zero
        
        ind=find([amplifier_channels.chip_channel]==i); % ind to find this chan
        
        if isempty(ind)
            % then did not recoprd this chan...
            disp(['PROBLEM - do not have data for neural chan (skipping) ' num2str(i)]);
            continue
        end
        
        chanstmp=[chanstmp i];
        dattmp=[dattmp; amplifier_data(ind, :)];
    end
    
    AmpChans_zero=chanstmp;
    amplifier_data=dattmp;
end



%% =========== PLOT ALL CHANNELS (RAW, no filter)
% --- plot initiate
figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

hsplots=[];


% ==== DIG CHANNELS




% ==== ANALOG (board aux) inputs
if isempty(AnalogChans_zero)
    
else
    for i=1:length(AnalogChans_zero)
        
        % == PLOT
        % 1) raw amplitude
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['analog chan: ' num2str(AnalogChans_zero(i))]);
        plot(tt, board_adc_data(i, :), 'b');
        hsplots=[hsplots hsplot];
        
        % 2) spectrogram
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['analog chan: ' num2str(AnalogChans_zero(i))]);
        [t, f, spec]=fn_extract_sound_spec(frequency_parameters, board_adc_data(i, :));
        imagesc(t, f, spec);
        axis([t(1) t(end) f(1) f(end)]);
        hsplots=[hsplots hsplot];
    end
end



% ==== AMPLIFIER CHANNELS (neural)
if isempty(AmpChans_zero)
    
else
    for i=1:length(AmpChans_zero)
        
        % == PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['neural chan: ' num2str(AmpChans_zero(i))]);
        plot(tt, amplifier_data(i, :), 'k');
        hsplots=[hsplots hsplot];
    end
end


% ========== SUBTRACTED NEURAL CHANNELS
% -- chan2 minus chan1
chan1=19; chan2=23;
ind1=AmpChans_zero==chan1;
ind2=AmpChans_zero==chan2;

dat=amplifier_data(ind2, :) - amplifier_data(ind1, :);

        % == PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['neural chan: ' num2str(chan2) ' minus ' num2str(chan1)]);
        plot(tt, dat, 'm');
        hsplots=[hsplots hsplot];




% -- link all
linkaxes(hsplots, 'x');



%% =========== PLOT ALL CHANNELS (filtering neural)
% --- plot initiate
figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

hsplots=[];

% ---- get timebins
if exist('amplifier_data', 'var')
    tt=1:length(amplifier_data(1, :));
    tt=tt./fs; % in sec
else
    tt=1:length(board_adc_data(1, :));
    tt=tt./fs; % in sec
end


% ==== DIG CHANNELS




% ==== ANALOG (board aux) inputs
if isempty(AnalogChans_zero)
    
else
    for i=1:length(AnalogChans_zero)
        
        % == PLOT
        % 1) raw amplitude
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['analog chan: ' num2str(i)]);
        plot(tt, board_adc_data(i, :), 'b');
        hsplots=[hsplots hsplot];
        
        % 2) spectrogram
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['analog chan: ' num2str(i)]);
        [t, f, spec]=fn_extract_sound_spec(frequency_parameters, board_adc_data(i, :));
        imagesc(t, f, spec);
        axis([t(1) t(end) f(1) f(end)]);
        hsplots=[hsplots hsplot];
    end
    
    
    
end



% ==== AMPLIFIER CHANNELS (neural)
if isempty(AmpChans_zero)
    
else
    for i=1:length(AmpChans_zero)
        
        % === filter
        dat=filtfilt(filt_b,filt_a,amplifier_data(i, :));
        
        % == PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['neural chan (filt): ' num2str(AmpChans_zero(i))]);
        plot(tt, dat, 'k');
        hsplots=[hsplots hsplot];
    end
    
    
end


% ==================== SUBTRACTED NEURAL CHANNELS
% -- chan2 minus chan1
chan1=19; chan2=23;
ind1=AmpChans_zero==chan1;
ind2=AmpChans_zero==chan2;

dat=amplifier_data(ind2, :) - amplifier_data(ind1, :);
        dat=filtfilt(filt_b,filt_a,amplifier_data(i, :)); % filter

        % == PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['neural chan: ' num2str(chan2) ' minus ' num2str(chan1)]);
        plot(tt, dat, 'm');
        hsplots=[hsplots hsplot];



% -- link all
linkaxes(hsplots, 'x');



end


function [t, f, sptemp] = fn_extract_sound_spec(frequency_parameters, rawdat)

fs=frequency_parameters.board_adc_sample_rate;
window = 0.016*fs; % make it 16ms, to match evtaf stuff
nfft=2^nextpow2(window); % go to next pow 2 size
olap=0.8;
noverlap=olap*window;


% ==== filter song
songdat=(rawdat);
filter_type='hanningfir';
F_low  = 500;
F_high = 8000;
songdat=bandpass(songdat,fs,F_low,F_high,filter_type);


% === colect spectrogram
[sp, f, t] = spectrogram(songdat, window, noverlap, nfft, fs);
t=t.*1000;
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


end





