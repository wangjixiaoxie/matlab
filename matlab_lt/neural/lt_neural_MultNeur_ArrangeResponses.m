function lt_neural_MultNeur_ArrangeResponses(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, WHOLEBOUTS_edgedur)

%% ======= for each neuron get mean response for one motif, and arrange by date, depth and channel

% motif_regexpr_str={'a(b)', 'j(b)'};
% motif_predur=0.1;
% motif_postdur=0.2;
% 
% LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
% WHOLEBOUTS_edgedur = '';

[MOTIFSTATS, MotifTime_med, spktimefield] = lt_neural_MultNeur_MotifRasters_v2(NeuronDatabase, motif_regexpr_str, ...
    motif_predur, motif_postdur, LinScaleGlobal, WHOLEBOUTS_edgedur);

%% smoothing window (neural)
window=0.02; windshift=0.004;
RandomSong =1; % set as default, for picking spectrogram

%% PLOT (for each channel, arrange with x = time; y = depth)

ChannelList = unique([NeuronDatabase.neurons.chan]);

for cc=1:length(ChannelList)
    ChanPlot = ChannelList(cc);
    lt_figure; hold on;
    
    % ------- list of "neurons"
    NeuronsToPlot = find([NeuronDatabase.neurons.chan] == ChanPlot);
    
    % ------ list of depths
    DepthsToPlot = [NeuronDatabase.neurons(NeuronsToPlot).electrode_depth];
    
    % ------- list of song times
    FirstSongTimeList_days = [];
    LastSongTimeList_days = [];
    FirstSongDatenumList = [];
    LastSongDatenumList = [];
    
    for nnn = NeuronsToPlot
       tmp = min([MOTIFSTATS.neurons(nnn).motif(1).SegmentsExtract.song_datenum]);
       FirstSongTimeList_days = [FirstSongTimeList_days tmp];        
       
       tmp = max([MOTIFSTATS.neurons(nnn).motif(1).SegmentsExtract.song_datenum]);
       LastSongTimeList_days = [LastSongTimeList_days tmp];
    end
    FirstSongDatenumList = FirstSongTimeList_days; %
    LastSongDatenumList = LastSongTimeList_days; %   
    
    firstday = datestr(min(FirstSongTimeList_days), 'ddmmmyyyy');
    tmp = lt_convert_EventTimes_to_RelTimes(firstday, FirstSongTimeList_days);
    earliesttime = min(tmp.FinalValue);
    FirstSongTimeList_days = tmp.FinalValue - earliesttime; %
    
    tmp = lt_convert_EventTimes_to_RelTimes(firstday, LastSongTimeList_days);
    LastSongTimeList_days = tmp.FinalValue - earliesttime; %
    
    % ---- prepare for subplots
    numcols = length(NeuronsToPlot);
    depthrange = max(DepthsToPlot) - min(DepthsToPlot);
    depthmin = min(DepthsToPlot);
    hsplots = [];
    
    for zzz = 1:length(NeuronsToPlot)
    % --- Prepare subplot
    ypos = 1-(1/numcols) - ((numcols-1)/numcols)*(DepthsToPlot(zzz) - depthmin)/depthrange;
    xpos = (zzz-1)*1/numcols;
    xwidth = 0.9*1/numcols;
    yheight = 0.9*1/numcols;
    
    multiplier = 0.85;
    ypos = 0.05+multiplier*ypos; xpos = 0.05+multiplier*xpos; xwidth = multiplier*xwidth; yheight = multiplier*yheight;
    
%     disp([xpos ypos xwidth yheight]);
    hsplot = lt_subplot(numcols, numcols, zzz);
    set(hsplot, 'Position', [xpos ypos xwidth yheight]); hold on;
    hsplots = [hsplots hsplot];
    
    title({[num2str(FirstSongTimeList_days(zzz)*24, '%3.2g') ' to ' num2str(LastSongTimeList_days(zzz)*24, '%3.2g')...
        ' hrs - Depth ' num2str(DepthsToPlot(zzz))],...
        [datestr(FirstSongDatenumList(zzz), 'ddmmm-HHMM') ' to ' ...
        datestr(LastSongDatenumList(zzz), 'ddmmm-HHMM')]});
    xlabel(['neuronID ' num2str(NeuronsToPlot(zzz))]);
    
    
    % ================================== PLOT RESPONSE HERE
    i = NeuronsToPlot(zzz);
    NumMotifs = length(MOTIFSTATS.neurons(i).motif);
    plotcols_motif = lt_make_plot_colors(NumMotifs, 0, 0);
    for m=1:NumMotifs
        if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
            % then there is no data for this neuron
            continue
        end
 
                segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end

        if length(Yspks)>1
        [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
        
        % --- plot
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols_motif{m}}, 1);
        end

                % -- line for motif onset
        line([motif_predur motif_predur], ylim, 'Color','k');
        if LinScaleGlobal==1
            line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','k'); % end
        end
        ylim([0 max(ymean_hz)+20]);

    end
    
    end
    axis tight
    linkaxes(hsplots, 'xy');
    
    for m=1:NumMotifs
    % ============= PLOT EXAMPLE SPECTROGRAM
    ypos = 1-(1/numcols) - ((numcols-1)/numcols)*(DepthsToPlot(1) - depthmin)/depthrange;
    ypos = ypos - 1.7*m*0.9*1/numcols;
    if ypos<0
    ypos = 1-(1/numcols) - ((numcols-1)/numcols)*(DepthsToPlot(1) - depthmin)/depthrange;
    ypos = ypos + 1.7*m*0.9*1/numcols;
    if ypos>1
        disp('PROBLEM !!! - not enough space to plot all motifs...');
    end
    end
    xpos = 0;
    xwidth = 0.9*1/numcols;
    yheight = 0.8*1/numcols;
    
    multiplier = 0.85;
    ypos = 0.05+multiplier*ypos; xpos = 0.05+multiplier*xpos; xwidth = multiplier*xwidth; yheight = multiplier*yheight;
    htmp = lt_subplot(numcols, numcols, 1);
    set(htmp, 'Position', [xpos ypos xwidth yheight]); hold on;

    songdat=[];
    fs=[];
    if LinScaleGlobal==1
        % --- find the trial that has close to median motif time
        for i=NeuronsToPlot
            if ~isfield(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            
            numtrials=length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract);
            
            for j=1:numtrials
                if abs(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur - ...
                        MotifTime_med)<0.01*MotifTime_med; % good enough, less than 1% of motif median.
                    
                    songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
                    fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
                    %                 motifdur=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).actualmotifdur;
                    break
                end
            end
        end
    end
    if RandomSong==1
        % make it empty so will pick random song
        songdat=[];
    end
    
    if isempty(songdat)
        % then pick a random song
        i=NeuronsToPlot(randi(length(NeuronsToPlot)));
        j=randi(length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract),1);
        songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).songdat;
        fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(j).fs;
    end
    
    % --- plot
    lt_plot_spectrogram(songdat, fs, 1, 0);
    line([motif_predur motif_predur], ylim, 'Color','w');
    if LinScaleGlobal==1
        line([motif_predur+MotifTime_med motif_predur+MotifTime_med], ylim, 'Color','w'); % end
    end
    if RandomSong==1
        title(['rand song: ' motif_regexpr_str{m}], 'Color', plotcols_motif{m});
    else
        title(['dur close to median: ' motif_regexpr_str{m}], 'Color', plotcols_motif{m});
    end
    
    end

    lt_subtitle(['CHANNEL ' num2str(ChanPlot)]);
    
end
