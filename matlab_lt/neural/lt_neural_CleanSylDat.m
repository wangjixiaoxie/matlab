function lt_neural_CleanSylDat(basedir, subdirs, Batchfiles, pretime, posttime, ...
    chanstoplot, skipifdone, plotspecgram, freqrange)
%% params

% ----- these three things, indices must be aligned.
% basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110517_RALMANlearn2/';
% subdirs = {'', ''}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
% Batchfiles = {'BatchSw1Pre', 'BatchSw1Post'}; % must be one for each subdir;

% pretime = 0.1; % sec, from onset
% posttime = 0.1; % sec, from offset
% chanstoplot = [9 14];
% skipifdone =1; % if 1, then skips if already done and old notmat labels are identical to current and chans to plot identical. if 0, then always redos
% plotspecgram=1; % then will plot spectrograms alongsize filtered neural.
% freqrange = [1800 4000]; % plots mean power over this range. if 0, then doesn't plot


%%
if plotspecgram==1
    % then will plot spectrograms alongsize filtered neural.
    lt_switch_chronux(1);
%     colormap('gray');
end
%% lt 11/27/17 - goes thru all syls for all songs and desired chans in batch, and asks you which to note down as noise
% --- will save a [filename].noise.mat file in the same directory. it
% contians syl positions for all syls that you call noisy.



% === go thru all song files. for each song file plot a figure, sorted by
% syllable, asking which files you want to remov



% =========== GO THRU ALL BATCHES
for i=1:length(subdirs)
    
    songdir = [basedir subdirs{i}];
    cd(songdir);
    % where the songs and batch are located
    batchf = Batchfiles{i}; % contains .rhd names
    
    % ================================
    fid = fopen(batchf);
    
    fline = fgetl(fid);
    while ischar(fline)
        
        % --- load stuff for this song
        if ~exist([fline '.not.mat'])
                        fline = fgetl(fid);
    continue
        end
        
        notdat = load([fline '.not.mat']);
        
        % -------------- check if already done
        if exist([fline '.noise.mat'], 'file');
            % --- check that old notmat matches
            noisetmp = load([fline '.noise.mat']);
            if all(noisetmp.trialsToRemove.notdat.labels == notdat.labels) ...
                    & all(find(~cellfun('isempty', noisetmp.trialsToRemove.chan)) == sort(chanstoplot)) ... % channels match
                    & skipifdone ==1
                disp(['SKIPPING ' fline ' -  previously done!!']);
                fline = fgetl(fid);
                continue
            end
        end
        
        [amplifier_data,~,frequency_parameters, ~, ...
            ~, amplifier_channels, ~, ~] =...
            pj_readIntanNoGui(fline);
        fs = frequency_parameters.amplifier_sample_rate;
        sylstoplot = unique(notdat.labels);
        
        % ------------- do one by one for each chan
        trialsToRemove = struct;
        trialsToRemove.chan = cell(1, 32);
        for chan = chanstoplot
            figcount=1;
            subplotrows=6;
            if ~isempty(freqrange) & plotspecgram==1
                subplotcols = 3;
            else
            subplotcols=2;
            end
            fignums_alreadyused=[];
            hfigs=[];
            
            PositionsCollected = [];
            count = 1;
            
            songdat = amplifier_data([amplifier_channels.chip_channel] == chan, :);
            
            for syl = sylstoplot
                
                if strcmp(syl,'-')
                    continue
                end
                
                listofpos = strfind(notdat.labels, syl);
                for pos=listofpos
                    
                    % ==================== extract dat and plot
                    % --- onset, offset
                    onset = notdat.onsets(pos);
                    onset = round(fs*onset/1000); % convert to samps
                    
                    onsetextract = onset -pretime*fs; % pretime and posttime
                    offsetextract = onset + posttime*fs;
                    
                    % --- extract neural
                    dat = songdat(onsetextract:offsetextract);
                    t = [1:length(dat)]./fs;
                    
                    % --- filter
                    datfilt = lt_neural_filter(dat, frequency_parameters);
                    
                    % --- plot
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['ch' num2str(chan) '-' syl ', #' num2str(count)]);
                    plot(t, datfilt, 'k');
                    axis tight
                    ylim([-200 200]);
                    line([pretime pretime], ylim);
                    
                    count = count+1;
                    
                    % --- get spectrogram
                    if plotspecgram ==1
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['ch' num2str(chan) '-' syl ', #' num2str(count)]);
                        
%                         specwin = 0.05; % sec
                        specwin = 0.02; % sec
                        specstep = 0.01;
                        params = struct;
%                         params.fpass = [0 1000/fs];
                        [S,t,f] = mtspecgramc(dat', [specwin specstep]*fs, params);
                        t = t./fs;
                        f = f.*fs;
                        %                     imagesc(t, f, S');
                        plot_matrix(S, t, f, 'l');
                        axis tight
                        ylim([0 4000]);
                        
                        % ----- average over freq range if desired
                        if ~isempty(freqrange)
                            
                            inds = f>=freqrange(1) & f<=freqrange(2);
                            Y = 10*log10(mean(S(:,inds),2));
                            
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['ch' num2str(chan) '-' syl ', #' num2str(count)]);
                        plot(t, Y, '-o');
                        axis tight
                        ylim([0 75]);
                        end
                                            
                    end
                end
                % --- append positions
                PositionsCollected = [PositionsCollected listofpos];
            end
            assert(length(PositionsCollected)== count-1, 'asdfasd');
            
            % ---- maximize all figures
            h = get(0,'children');
            for kk = 1:length(h)
                set(h(kk),'units','normalized','outerposition',[0 0 1 1])
                %             set(h(kk), 'Position', [0.0217 0 0.338 0.959]);
            end
            
            tt = input([fline ', Which trials to remove (e.g. 1 2 3 or 5:10)? '], 's');
            tt = str2num(tt);
            disp(['(to remove) trials: ' num2str(tt)]);
            trialsToRemove.chan{chan} = sort(PositionsCollected(tt)); % convert from trials to actual position in labels
            if isempty(tt)
                trialsToRemove.chan{chan} = nan;
            end
            close all;
            
        end
        % ====================== save trials to remove for this song file
        fnamesave = [fline '.noise.mat'];
        trialsToRemove.dateCleaned = lt_get_timestamp(0);
        trialsToRemove.pretime = pretime;
        trialsToRemove.posttime = posttime;
        trialsToRemove.notdat = notdat;
        save(fnamesave, 'trialsToRemove');
        disp('SAVED!');
        
        % -- load next file
        fline = fgetl(fid);
        
    end
end

fclose('all');
%%

lt_switch_chronux(0);

%% sanity check, plot stuff that was saved - verify that all are noisy.
if (0)
    
    lt_figure; hold off;
    
    chan = 14;
    fname = 'pu69wh78_171105_111100.rhd';
    [amplifier_data,~,frequency_parameters, ~, ...
        ~, amplifier_channels, ~, ~] =...
        pj_readIntanNoGui(fname);
    notdat = load([fname '.not.mat']);
    
    tmp = load([fname '.noise.mat']);
    % tmp.trialsToRemove.chan{chan} = 44:48
    for i=1:length(tmp.trialsToRemove.chan{chan})
        pos = tmp.trialsToRemove.chan{chan}(i);
        ons = fs*notdat.onsets(pos)/1000;
        
        dat = amplifier_data([amplifier_channels.chip_channel]==chan, :);
        dat = dat(ons-0.1*fs:ons+0.1*fs);
        dat = lt_neural_filter(dat, frequency_parameters);
        
        plot(dat);
        
        pause
    end
end