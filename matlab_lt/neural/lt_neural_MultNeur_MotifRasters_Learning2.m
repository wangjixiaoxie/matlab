function lt_neural_MultNeur_MotifRasters_Learning2(NeuronDatabase, motif_regexpr_str, motif_predur, ...
    motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, TrialBinSize, DivisorBaseSongs)
%% constrain to 1 neuron, as WNchangeDateStrings is for only one expt - and lots of figures


if length(NeuronDatabase.neurons)>1
    NeuronDatabase.neurons=NeuronDatabase.neurons(1);
end

% -- 1st ind is always when WN was turned on. other inds is for switches,
% etc.
WNchangeDateStrings{1}=NeuronDatabase.neurons(1).LEARN_WNonDatestr;
WNchangeDateStrings=[WNchangeDateStrings NeuronDatabase.neurons(1).LEARN_WNotherImportantDates];


%% for each motif and each neuron, make one plot, showing progression over time of
% motif_regexpr_str={'g(h)'};
% motif_predur=0.2;
% motif_postdur=0.1;
% LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
%
%
% FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
% FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
%     % +1 is 1 after token
% FFparams.FF_sylName='h'; % Optional: what syl do you expect this to be? if incompatible will raise error
%     % not required (can set as []);
% FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
%             'v', [2450 4300]};
% % FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
% %             'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
% FFparams.cell_of_FFtimebins={'h', [0.034 0.038], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % WN on g H
% % NOTE: will also determine whether was hit or miss, based on WN sound
% % detection.
%
% % LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'};
%
% OnlyPlotNoHit=1; % then only plots trials that were not hit (WN)




%% PARAMS
NumMotifs=length(motif_regexpr_str);
NumNeurons=length(NeuronDatabase.neurons);

% smoothing window (neural)
window=0.01; windshift=0.002;

plotTimesOnRaster=1; % hh:mm on raster plot

determineTrialBinBasedOnBaselineN=1; % if 1, then chooses binsize based on
% baseline num songs (i.e. baseNumSongs/DivisorBaseSongs); if 0, then = TrialBinSize
% note: assumes WNchangeDateStrings{1} is transition from baseline to WN
% on.
% DivisorBaseSongs=1; % divide base N with this
% TrialBinSize=10;


%% EXTRACT DATA (for each neuron x motif)
MOTIFSTATS=struct;
for i=1:NumNeurons
    try
    cd(NeuronDatabase.global.basedir);
    catch err
    cd(NeuronDatabase.neurons(i).basedir);    
    end
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(i).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);
    
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    extractSound = 1;
    tic
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractSound);
    toc
    
    % --- EXTRACT DAT
    % - do this one time for each desired motif
    for j=1:NumMotifs
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        keepRawSongDat=1;
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, keepRawSongDat);
        
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
    end
end

% doclose=input('type "c" to close figures (otherwise ENTER)', 's');
doclose='c';
if strcmp(doclose, 'c');
    close all;
end
%
% disp('PRESS ANYTHING TO CONTINUE');
% pause;
% close all;
%
%% === MAKE SURE ALL TRIALS ARE IN TEMPORAL ORDER

for i=1:NumNeurons
    for m=1:NumMotifs
        
        segextract = MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        
        all_datenums=[segextract.song_datenum];
        
        [~, inds] = sort(all_datenums);
        
        if all(diff(inds))~=1
            disp('HAD TO REORDER !! ---- not a problem');
            
            MOTIFSTATS.neurons(i).motif(m).SegmentsExtract=segextract(inds);
        end
    end
end

%% ===== PLOT, FOR EACH MOTIF AND EACH NEURON, PLOT SUMMARY PLOT
PlotSmallTimeWindow = 1; % then cuts off end of motif. keeps only based on what is specified in postdur.

spktimefield='spk_Times';
plotcols=lt_make_plot_colors(NumNeurons, 0, 0);

for i=1:NumNeurons
    
    for m=1:NumMotifs
        
        lt_figure; hold on;
        hsplots=[];
        hsplots2=[];
        
        numrows = 10;
        
        % ===== SONG SPECTROGRAM
        hsplot= lt_subplot(numrows, 5, 2:3); hold on; hsplots =[hsplots hsplot];
        % then pick a random song
        tmp=randi(length(MOTIFSTATS.neurons(i).motif(m).SegmentsExtract),1);
        songdat=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(tmp).songdat;
        fs=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract(tmp).fs;
        % --- plot
        lt_plot_spectrogram(songdat, fs, 1, 0);
        %        line([motif_predur motif_predur], ylim, 'Color','w');
        
        axis tight
        if PlotSmallTimeWindow==1
            xlim([0 motif_predur+0.15]);
        end
        
        
        % ===== RASTERS
        subplotinds=sort([7:5:50 8:5:50]);
        hsplot = lt_subplot(numrows, 5, subplotinds); hold on; hsplots =[hsplots hsplot];
        yind=1; % will increment to plot all
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        params=MOTIFSTATS.neurons(i).motif(m).Params;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        % --- IF ONLY PLOT NO HIT, THEN FIRST PULL OUT THOSE TRIALS
        if OnlyPlotNoHit==1
            indsWithNoHits=[segextract.hit_WN]~=1;
            segextract=segextract(indsWithNoHits);
        end
        
        % --- PLOT DAT
        numtrials=length(segextract);
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            
            if PlotSmallTimeWindow ==1
               windowtmp = [0 motif_predur+0.15];
                spktimes = spktimes(spktimes>windowtmp(1) & spktimes<windowtmp(2));
            end
            
            if ~isempty(spktimes)
                for jj=1:length(spktimes)
                    spk=spktimes(jj);
                    line([spk spk], -[j-0.4 j+0.4], 'Color', plotcols{i}, 'LineWidth', 1);
                end
            end
            
            % ==== plot time of trial
            if plotTimesOnRaster==1
                trialtime=segextract(j).song_datenum;
                trialtime=datestr(trialtime, 'HH:MM');
                lt_plot_text(0, -j, trialtime, [0.4 0.4 0.4]);
            end
            
            % ========= WAS TRIAL A HIT?
            if segextract(j).hit_WN==1
                %                 line([-0.02 -0.002], [-j -j], 'Color','k', 'LineWidth' , 2);
                % note when WN hit
                X=[segextract(j).WNonset_sec  segextract(j).WNoffset_sec ...
                    segextract(j).WNoffset_sec  segextract(j).WNonset_sec];
                Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
                h=patch(X, Y, 'w', 'FaceColor', [0 0 0], 'EdgeColor', [0.8 0.8 0.8]);
                set(h, 'FaceAlpha', 0.1)
                
            end
            
        end
        line([motif_predur motif_predur], ylim)
        axis tight
        
        % ====== ANNOTATE WN CHANGE TIMEPOINTS
        for dd=1:length(WNchangeDateStrings)
            datestring=WNchangeDateStrings{dd};
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            
            allsongdnums=[segextract.song_datenum];
            sylind=find(allsongdnums>dnum, 1, 'first');
            
            line(xlim, [-sylind+0.5 -sylind+0.5], 'Color', 'k');
            lt_plot_text(0, -sylind+0.5, datestring, 'k');
        end
        
        
%         % ===== RUNNING AVERAGE FOR ENTIRE DAY
%         hsplot=lt_subplot(8,4,29:31); hold on; hsplots=[hsplots hsplot];
%         
%         numtrials=length(segextract);
%         Yspks={};
%         for j=1:numtrials
%             inds=segextract(j).spk_Clust==clustnum;
%             spktimes=segextract(j).(spktimefield)(inds);
%             Yspks{j}=spktimes;
%         end
%         % -- convert to smoothed rate
%         [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
%         % --- plot
%         shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
%         lt_plot_zeroline
%         ylim([0 300]);
        
        
        % ===== RUNNING AVERAGE OF SMOOTHED FIRING RATE (OVER HANDFUL OF TRIALS
        % --- average over every N trials
        if determineTrialBinBasedOnBaselineN==1
            datestring=WNchangeDateStrings{1}; % assumes this is WN on.
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            N=ceil(sum([segextract.song_datenum] < dnum)/DivisorBaseSongs); % num, baseline rends divid by divisor.
        else
            N=15;
        end
        
        for j=1:ceil(numtrials/N);
            trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
            
            hsplot=lt_subplot(ceil(numtrials/(N)), 5, [5*(j)-1 5*j]); hold on; hsplots=[hsplots hsplot];
            hsplots2=[hsplots2 hsplot];
            lt_plot_text(0, 250, ['trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'b');
            Yspks={};
            for jj=trialsToPlot
                inds=segextract(jj).spk_Clust==clustnum;
                spktimes=segextract(jj).(spktimefield)(inds);
                
                            if PlotSmallTimeWindow ==1
               windowtmp = [0 motif_predur+0.15];
                spktimes = spktimes(spktimes>windowtmp(1) & spktimes<windowtmp(2));
                            end

            
                Yspks=[Yspks spktimes];
            end
            % -- convert to smoothed rate
            [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            
            % -- save first bin to overlay over other bins
            if j==1;
                xbin1=xbin;
                ymean1=ymean_hz;
                ysem1=ysem_hz;
            else
                shadedErrorBar(xbin1, ymean1, ysem1, {'Color', [0.75 0.75 0.75]}, 1);
            end
            % --- plot
            if length(trialsToPlot)>1 % mean
                shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
                lt_plot_zeroline
                ylim([0 300]);
                line([motif_predur motif_predur], ylim);
            end
        end
        
        
        % === PLOT FF
        % TRAJECTORY
        FFvals=[segextract.FF_val];
        Tvals=[segextract.song_datenum];
        if ~isempty(FFvals)
            subplotinds=sort([6:5:50]);
            lt_subplot(numrows,5,subplotinds); hold on; % by rendition
            
            ylabel('rendition');
            xlabel('FF');
            plot(FFvals, -[1:length(FFvals)], 'o', 'Color', plotcols{i});
            
            % WN change points
            for dd=1:length(WNchangeDateStrings)
                datestring=WNchangeDateStrings{dd};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                allsongdnums=[segextract.song_datenum];
                sylind=find(allsongdnums>dnum, 1, 'first');
                
                line(xlim, -[sylind-0.5 sylind-0.5], 'Color', 'k');
                lt_plot_text(FFvals(1), -(sylind+0.5), datestring, 'k');
            end
            
            % ----- overlay +/- 1 SD of baseline
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                allsongdnums=[segextract.song_datenum];
                baseinds=allsongdnums<dnum;
                 baseFF_mean = mean(FFvals(baseinds));
                 baseFF_STD = std(FFvals(baseinds));
                 
                 line([baseFF_mean baseFF_mean], ylim, 'Color', 'k');
                 line([baseFF_mean-baseFF_STD baseFF_mean-baseFF_STD], ...
                     ylim, 'Color', 'k', 'LineStyle', '--');
                 line([baseFF_mean+baseFF_STD baseFF_mean+baseFF_STD], ...
                     ylim, 'Color', 'k', 'LineStyle', '--');
                 
% ==== overlay mean FF for the windows used to extract FR
        RendMean = [];
        FFmeans = [];
        FFsems = [];
        FFstds = [];
        for j=1:ceil(numtrials/N);
            trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
            
            ffvalstmp = FFvals(trialsToPlot);
            
            ffmean = mean(ffvalstmp);
            ffsem = lt_sem(ffvalstmp);
            ffstd = std(ffvalstmp);
            
            RendMean = [RendMean mean(trialsToPlot)];
            FFmeans = [FFmeans ffmean];
            FFsems = [FFsems ffsem];
            FFstds = [FFstds ffstd];
         end
         plot(FFmeans, -RendMean, 'o-', 'MarkerSize',8, 'Color', 'k', 'LineWidth', 2);  
         for j=1:length(FFstds)
            line([FFmeans(j)-FFstds(j) FFmeans(j)+FFstds(j)], ...
                -[RendMean(j) RendMean(j)], 'Color', 'k', 'LineWidth', 2);
         end
            
         axis tight
%             lt_subplot(1,2,2); hold on; % by time
%             xlabel('time');
%             plot(Tvals, FFvals, 'o', 'Color', plotcols{i});
%             datetick('x', 'ddmmm-HHMM', 'keepticks');
%             % WN change points
%             for dd=1:length(WNchangeDateStrings)
%                 datestring=WNchangeDateStrings{dd};
%                 dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
%                 
%                 line([dnum dnum], ylim, 'Color', 'k');
%                 %                 lt_plot_text(sylind+0.5, FFvals(1), datestring, 'k');
%             end
%             lt_subtitle(motif_regexpr_str{m})
            
        end
        
                % === all subplots
        linkaxes(hsplots, 'x');
        linkaxes(hsplots2, 'xy');
        lt_subtitle(['neuron ' num2str(i) ', motif: ' motif_regexpr_str{m}]);

        
    end
    
end

%% ONE PLOT FOR EACH NEURON (SHOW ALL MOTIFS)








%% ===== SUMMARY PLOT, PLOT RUNNING AVERAGE OF Z-SCORE RELATIVE TO BASELINE
% PLOT VS. TIME OR VS. REND

if (0)
%% ===== PLOT, FOR EACH MOTIF AND EACH NEURON, PLOT SUMMARY PLOT
spktimefield='spk_Times';
plotcols=lt_make_plot_colors(NumNeurons, 0, 0);

for i=1:NumNeurons
    
    for m=1:NumMotifs
        
        
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        % --- IF ONLY PLOT NO HIT, THEN FIRST PULL OUT THOSE TRIALS
        if OnlyPlotNoHit==1
            indsWithNoHits=[segextract.hit_WN]~=1;
            segextract=segextract(indsWithNoHits);
        end
        
        %         % ====== ANNOTATE WN CHANGE TIMEPOINTS
        %         for dd=1:length(WNchangeDateStrings)
        %             datestring=WNchangeDateStrings{dd};
        %             dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
        %
        %             allsongdnums=[segextract.song_datenum];
        %             sylind=find(allsongdnums>dnum, 1, 'first');
        %
        %             line(xlim, [-sylind+0.5 -sylind+0.5], 'Color', 'k');
        %             lt_plot_text(0, -sylind+0.5, datestring, 'k');
        %         end
        %
        
        % ===================== collect spikes for all trials
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        % ===================== BASELINE mean and std spiking
        datestring=WNchangeDateStrings{1}; % assumes this is WN on.
        dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
        inds=[segextract.song_datenum] < dnum; %
        [xbase, ~, ~, ymean_hz_base, ~, ~, ystd_hz_base] = ...
            lt_neural_plotRastMean(Yspks(inds), window, windshift, 0, '');
        
        
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % ===================== OVERLAPPING SMOOTHED d-prime rate
        lt_figure; hold on;
        title(['neuron' num2str(i) ': ' motif_regexpr_str{m}]);
        xlabel('time (within trial');
        ylabel('bin median time');
        zlabel('d-prime vs. base');
        
        cmap=cool((numtrials-TrialBinSize+1));
        useRealTime=1;
        AllDprimeTraj={};
        AllTrialsPreBase=[];
        AllTrialsDurWN=[];
        AllMedianTimes=[];
        
        for j=1:(numtrials-TrialBinSize+1);
            trialbins=j:(j+TrialBinSize-1);
            
            %  -- figure out if all trials are dur baseline or dur WN or
            %  micture.
            datestring=WNchangeDateStrings{1}; % assumes this is WN on.
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            if segextract(trialbins(end)).song_datenum < dnum
                AllTrialsPreBase=[AllTrialsPreBase 1];
                AllTrialsDurWN=[AllTrialsDurWN 0];
            elseif segextract(trialbins(1)).song_datenum >= dnum
                AllTrialsPreBase=[AllTrialsPreBase 0];
                AllTrialsDurWN=[AllTrialsDurWN 1];
            else
                AllTrialsPreBase=[AllTrialsPreBase 0];
                AllTrialsDurWN=[AllTrialsDurWN 0];
            end
            
            [xtime, ~, ~, ymean_hz, ~, ~, ystd_hz] = lt_neural_plotRastMean(Yspks(trialbins), window, windshift, 0, '');
            
            % -- calc d-prime
            xinds=1:min([length(ymean_hz) length(ymean_hz_base)]);
            numer_tmp=ymean_hz(xinds)-ymean_hz_base(xinds);
            denom_tmp=sqrt(0.5*(ystd_hz_base(xinds).^2 + ystd_hz(xinds).^2));
            dprime=numer_tmp./denom_tmp;
            
            % -- get median time for this bin
            median_datenum=median([segextract(trialbins).song_datenum]);
            
            % days from start of expt
            firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
            eventtime=datestr(median_datenum, 'ddmmmyyyy-HHMM');
            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime});
            
            % ====== ADD LINE TO 3D PLOT
            Z=dprime;
            if useRealTime==1
                Y=ones(1, length(dprime)) * tmp.FinalValue;
            else
                Y=ones(1, length(dprime)) * j;
            end
            X=xtime(xinds);
            %             lt_plot_stem3(X, Y, Z, cmap(j, :), 1);
            plot3(X, Y, Z, 'Color', cmap(j, :));
            
            % --- add plane for 0, -1, 1
            
            % ====== collect all d-prime trajectories for heat map plot
            AllDprimeTraj=[AllDprimeTraj dprime];
            AllMedianTimes=[AllMedianTimes tmp.FinalValue];
        end
        
        axis tight
        zlim([-3 3]);
        set(gca, 'YDir', 'reverse')
        if (1) % plots grid mesh
            x=xlim; x=linspace(x(1), x(2), 10);
            y=ylim; y=linspace(y(1), y(2), 10);
            z=1*ones(length(x), length(y));
            hmesh=mesh(x, y, z);
            
            set(hmesh, 'FaceAlpha', 0.5)
            set(hmesh, 'EdgeColor', 'r')
            set(hmesh, 'LineWidth', 0.3)
            
            x=xlim; x=linspace(x(1), x(2), 10);
            y=ylim; y=linspace(y(1), y(2), 10);
            z=-1*ones(length(x), length(y));
            hmesh=mesh(x, y, z);
            
            set(hmesh, 'FaceAlpha', 0.5)
            set(hmesh, 'EdgeColor', 'g')
            set(hmesh, 'LineWidth', 0.3)
        end
        % -- plot line for when WN start
        if useRealTime==1
            datestring=WNchangeDateStrings{1}; % assumes this is WN on.
            firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {datestring});
            line(xlim, [tmp.FinalValue tmp.FinalValue], zlim, 'Color','k');
        end
        
        
        % =========== ALTERNATIVELY, PLOT HEAT MAP ACROSS TIME
        lt_figure; hold on;
        title(['neuron' num2str(i) ': ' motif_regexpr_str{m}]);
        xmax=min(cellfun(@length, AllDprimeTraj));
        AllDprimeMat=nan(length(AllDprimeTraj), xmax);
        for j=1:length(AllDprimeTraj)
            AllDprimeMat(j, :) = AllDprimeTraj{j}(1:xmax);
        end
        clear AllDprimeTraj
        colormap('spring');
        colorbar;
        imagesc(xtime(1:xmax), 1:size(AllDprimeMat,1), AllDprimeMat);
        
        % line for bins pre NW and post WN
        line(xlim, [sum(AllTrialsPreBase)+0.5 ...
            sum(AllTrialsPreBase)+0.5], 'Color','k');
        line(xlim, [sum(~AllTrialsDurWN)-0.5 ...
            sum(~AllTrialsDurWN)-0.5], 'Color','k');
        
        line(xlim, [TrialBinSize+0.5 TrialBinSize+0.5], 'Color','k','LineStyle', '--');
        lt_plot_text(0.01, TrialBinSize+0.6, 'nonoverlap with base bin', [0.4 0.4 0.4]);
        % ==== plot time of trial
        if plotTimesOnRaster==1
            for j=1:length(AllMedianTimes)
                lt_plot_text(0, j, num2str(AllMedianTimes(j),'%3.3g'), [0.4 0.4 0.4]);
            end
        end
        
        
    end
    
    
    %% ==== save figures
    cd(NeuronDatabase.global.basedir);

    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(i).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);
    
    try cd('FIGS')
    catch err
        mkdir('FIGS')
        cd('FIGS')
    end
    
    tstamp = lt_get_timestamp(0);
    
    mkdir(['MotifRasters_Learning_' tstamp]);
    cd(['MotifRasters_Learning_' tstamp]);
    lt_save_all_figs;    
    
end
end

%% === SUMMARIZE BY TAKING MEAN DPRIME OVER PREMOTOR WINODW
if (0) % SKIP, SINCE DOES THIS IN SUMMARY CODE.
    for i=1:NumNeurons
        figcount=1;
        subplotrows=6;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots1=[];
        hsplots2=[];
        
        
        for m=1:NumMotifs
            
            
            segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
            clustnum=NeuronDatabase.neurons(i).clustnum;
            
            % --- IF ONLY PLOT NO HIT, THEN FIRST PULL OUT THOSE TRIALS
            if OnlyPlotNoHit==1
                indsWithNoHits=[segextract.hit_WN]~=1;
                segextract=segextract(indsWithNoHits);
            end
            
            
            % ===================== collect spikes for all trials
            numtrials=length(segextract);
            Yspks={};
            for j=1:numtrials
                inds=segextract(j).spk_Clust==clustnum;
                spktimes=segextract(j).(spktimefield)(inds);
                Yspks{j}=spktimes;
            end
            
            % ===================== BASELINE mean and std spiking
            datestring=WNchangeDateStrings{1}; % assumes this is WN on.
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            inds=[segextract.song_datenum] < dnum; %
            [xbase, ~, ~, ymean_hz_base, ~, ~, ystd_hz_base] = ...
                lt_neural_plotRastMean(Yspks(inds), window, windshift, 0, '');
            
            
            % ========================= GET RUNNING D-PRIME OVER TIME, IN
            % SPECIFIC PREMOTOR WINDOW
            premotor_wind=[-0.08 0.02]; % in sec, relative to onset of token in motif.
            DprimeVals=[];
            AbsDprimeVals=[];
            MedianDatenumVals = [];
            MedianDayVals=[];
            FFVals=[];
            AreAllTrialsBaseline = [];
            AreAllTrialsDurWN = [];
            
            
            for j=1:(numtrials-TrialBinSize+1);
                trialbins=j:(j+TrialBinSize-1);
                
                [x, ~, ~, ymean_hz, ~, ~, ystd_hz] =  ...
                    lt_neural_plotRastMean(Yspks(trialbins), window, windshift, 0, '');
                
                % -- calc d-prime
                xinds=1:min([length(ymean_hz) length(ymean_hz_base)]);
                numer_tmp=ymean_hz(xinds)-ymean_hz_base(xinds);
                denom_tmp=sqrt(0.5*(ystd_hz_base(xinds).^2 + ystd_hz(xinds).^2));
                dprime=numer_tmp./denom_tmp;
                
                % -- slice dprime within premotor window
                xinds_to_keep = x>=(motif_predur + premotor_wind(1) + window/2) ...
                    & x<=(motif_predur + premotor_wind(2) - window/2); % accounting for smoothing window size.
                
                dprime_val = mean(dprime(xinds_to_keep));
                DprimeVals=[DprimeVals dprime_val];
                
                abs_dprime_val = mean(abs(dprime(xinds_to_keep)));
                AbsDprimeVals=[AbsDprimeVals abs_dprime_val];
                
                % -- get median time for this bin
                median_datenum=median([segextract(trialbins).song_datenum]);
                firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
                eventtime=datestr(median_datenum, 'ddmmmyyyy-HHMM');
                tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime}); % days from start of expt
                
                MedianDatenumVals = [MedianDatenumVals median_datenum];
                MedianDayVals = [MedianDayVals tmp.FinalValue];
                
                % -- FF
                mean_FF=mean([segextract(trialbins).FF_val]);
                FFVals=[ FFVals mean_FF];
                
                % -- all trials baseline?
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                if segextract(trialbins(end)).song_datenum < dnum
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 1]; % all trials base
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
                elseif segextract(trialbins(1)).song_datenum >= dnum
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 0];
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 1];
                else
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 0];
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
                end
                
                
            end
            
            
            % ==== PLOT
            if (0) % set as 1 to plot mean dprime, otherwise plots mean abs(dprime)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('mean premotor dprime');
                plot(MedianDayVals, DprimeVals, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), DprimeVals(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), DprimeVals(logical(AreAllTrialsDurWN)), 'or');
            else
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('mean abs(premotor dprime)');
                plot(MedianDayVals, AbsDprimeVals, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), AbsDprimeVals(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), AbsDprimeVals(logical(AreAllTrialsDurWN)), 'or');
            end
            
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots2=[hsplots2 hsplot];
            title(motif_regexpr_str{m}); ylabel('ff (hz)');
            plot(MedianDayVals, FFVals, 'sk');
            plot(MedianDayVals(logical(AreAllTrialsBaseline)), FFVals(logical(AreAllTrialsBaseline)), 'sb');
            plot(MedianDayVals(logical(AreAllTrialsDurWN)), FFVals(logical(AreAllTrialsDurWN)), 'sr');
            
            
        end
        
        linkaxes(hsplots1, 'xy');
        linkaxes(hsplots2, 'xy');
        lt_subtitle(['neuron' num2str(i)]);
    end
end

