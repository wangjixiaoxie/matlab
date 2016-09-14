function lt_neural_concatExplore_v2(batchf, ChansToPlot)

%% lt 8/17/16 - exploratory only - concats neural (all chans) + song
% does not save anything
%
% ChansToPlot=neural chans (chip chan)

    windowsize=0.03; % from -2sd to +2sd [for smoothed rectified]

%% put names of all files in batch into a cell array

filenames=lt_batchsong_NamesToCell(batchf);


%%  concatenate channels as vectors - save metadata information along the way

fs_all=[];
ampDat_all=cell(length(ChansToPlot),1); % chan x data
songDat_all=[];
transsamps=[];

% ----- COLLECT DATA
for i=1:length(filenames)
    
    fname=filenames{i};
    
    [amplifier_data,~ ,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(fname);
    
    % -- collect song
    songDat_all=[songDat_all board_adc_data(1,:)];
    transsamps=[transsamps length(board_adc_data(1,:))];
    
    % --- collect all amp dat
    for j=1:length(ChansToPlot)
        chan=ChansToPlot(j);
        ind=[amplifier_channels.chip_channel]==chan;
        dattmp=amplifier_data(ind, :); % dat for this chan
        
        ampDat_all{j}=[ampDat_all{j} dattmp]; % concat
        
        % --- collect all fs
        fs_all=[fs_all frequency_parameters.amplifier_sample_rate];
    end
    
    
end


% --- confirm that all have same sample rate
if length(fs_all)>1
    assert(all(diff(fs_all)==0), 'problem, FS not same for all file')
end

%% PLOT
figcount=1;
subplotrows=2*length(ChansToPlot)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];

% - song, raw
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots=[hsplots hsplot];
tt=[1:length(songDat_all)]/fs_all(1);
plot(tt, songDat_all);

%  line for each transition
transsamps=cumsum(transsamps);
for i=1:length(transsamps)
    x=tt(transsamps(i));
    line([x x], ylim)
end

% - plot each chan
for i=1:length(ChansToPlot)
    chan=ChansToPlot(i);
    
    % === 1) neural, filtered
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    title(['chan ' num2str(chan)]);
    
    datfilt=lt_neural_filter(ampDat_all{i}, frequency_parameters);
    plot(tt, datfilt, 'k');
    
    % === 2) neural, smoothed
    % === rectify
    datfilt=abs(datfilt);
    
    % == smoooth
    % Construct a gaussian window
    sigma=(windowsize/4)*fs_all(1); %
    numsamps=4*sigma; % (get 2 std on each side)
    alpha= numsamps/(2*sigma); % N/2sigma
    gaussFilter = gausswin(numsamps, alpha);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    dat_smrect = conv(datfilt, gaussFilter);
    
    % == PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    plot(tt, dat_smrect(numsamps/2:end-numsamps/2), 'k');
    hsplots=[hsplots hsplot];
    
    
    
    
end




linkaxes(hsplots, 'x');

