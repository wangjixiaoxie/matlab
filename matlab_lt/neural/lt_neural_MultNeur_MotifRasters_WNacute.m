function lt_neural_MultNeur_MotifRasters_WNacute(NeuronDatabase, motif_regexpr_str, ...
    motif_predur,motif_postdur, FFparams, NumTrialsBin, LinScaleGlobal)


%% PARAMS
NumMotifs=length(motif_regexpr_str);
NumNeurons=length(NeuronDatabase.neurons);

% smoothing window (neural)
window=0.02; windshift=0.004;

plotTimesOnRaster=1; % hh:mm on raster plot

determineTrialBinBasedOnBaselineN=1; % if 1, then chooses binsize based on
% baseline num songs (i.e. baseNumSongs/DivisorBaseSongs); if 0, then = TrialBinSize
% note: assumes WNchangeDateStrings{1} is transition from baseline to WN
% on.
DivisorBaseSongs=2;
TrialBinSize=10;

% NumTrialsBin=15;


%% EXTRACT DATA (for each neuron x motif)
MOTIFSTATS=struct;
for i=1:NumNeurons
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    % --- EXTRACT DAT
    % - do this one time for each desired motif
    for j=1:NumMotifs
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams);
        
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
    end
end

% doclose=input('type "c" to close figures (otherwise ENTER)', 's');
% if strcmp(doclose, 'c');
%     close all;
% end
%
% disp('PRESS ANYTHING TO CONTINUE');
% pause;
% close all;

%% ====== LINEARLY STRETCH IF NEEDED
suppresssomeplots=1;
if LinScaleGlobal==2
    
    spktimefield='spk_Times_scaled';
    for i=1:NumNeurons
        for j=1:NumMotifs
            if ~isfield(MOTIFSTATS.neurons(i).motif(j).SegmentsExtract, 'spk_Times')
                % then there is no data for this neuron
                continue
            end
            segmentdat=MOTIFSTATS.neurons(i).motif(j).SegmentsExtract;
            prms=MOTIFSTATS.neurons(i).motif(j).Params;
            [segmentdatout, prmsout] = lt_neural_LinTimeWarp(segmentdat, ...
                prms, []);
            
            if suppresssomeplots==1
                close gcf
                close gcf % close 2 figs
            end
            MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=segmentdatout;
            MOTIFSTATS.neurons(i).motif(j).Params=prmsout;
        end
    end
    
elseif LinScaleGlobal==0
    % THEN NO SCALING AT ALLL!
    spktimefield='spk_Times';
end

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


%% ===== PLOT, FOR EACH MOTIF AND EACH NEURON, PLOT RASTER (SEPARATING WN AND NO WN TRIALS
% spktimefield='spk_Times';
plotcols=lt_make_plot_colors(NumNeurons, 0, 0);

for i=1:NumNeurons
    
    clear WNchangeDateStrings;
    WNchangeDateStrings{1}=NeuronDatabase.neurons(i).LEARN_WNonDatestr;
    WNchangeDateStrings=[WNchangeDateStrings NeuronDatabase.neurons(i).LEARN_WNotherImportantDates];
    
    
    for m=1:NumMotifs
        
        hfig1= lt_figure; hold on;
        hfig2=lt_figure; hold on;
        figure(hfig1);
        hsplots=[];
        hsplots2=[];
        
        clustnum=NeuronDatabase.neurons(i).clustnum;
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        
        
        % ========================================================================= WN HIT TRIALS
        segextract_tmp=segextract([segextract.hit_WN]==1);
        
        
        % ===== SONG SPECTROGRAM
        hsplot= lt_subplot(8, 4, 1:2); hold on; hsplots =[hsplots hsplot];
        % then pick a random song
        tmp=randi(length(segextract_tmp),1);
        songdat=segextract_tmp(tmp).songdat;
        fs=segextract_tmp(tmp).fs;
        % --- plot
        lt_plot_spectrogram(songdat, fs, 1, 0);
        %        line([motif_predur motif_predur], ylim, 'Color','w');
        axis tight
        
        % ====== RASTER
        subplotinds=sort([5:4:21 6:4:22]);
        hsplot = lt_subplot(8, 4, subplotinds); hold on; hsplots =[hsplots hsplot];
        hsplots2=[hsplots2 hsplot];
        %         yind=1; % will increment to plot all
        
        numtrials=length(segextract_tmp);
        minsongsamps=[];
        for j=1:numtrials
            inds=segextract_tmp(j).spk_Clust==clustnum;
            spktimes=segextract_tmp(j).(spktimefield)(inds);
            
            if ~isempty(spktimes)
                for jj=1:length(spktimes)
                    spk=spktimes(jj);
                    line([spk spk], -[j-0.4 j+0.4], 'Color', plotcols{i}, 'LineWidth', 1);
                end
            end
            
            % ==== plot time of trial
            if plotTimesOnRaster==1
                trialtime=segextract_tmp(j).song_datenum;
                trialtime=datestr(trialtime, 'HH:MM');
                lt_plot_text(0, -j, trialtime, [0.4 0.4 0.4]);
            end
            
            % ========= WAS TRIAL A HIT?
            if segextract_tmp(j).hit_WN==1
                %                 line([-0.02 -0.002], [-j -j], 'Color','k', 'LineWidth' , 2);
                % note when WN hit
                X=[segextract_tmp(j).WNonset_sec  segextract_tmp(j).WNoffset_sec ...
                    segextract_tmp(j).WNoffset_sec  segextract_tmp(j).WNonset_sec];
                Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
                h=patch(X, Y, 'w', 'FaceColor', [0 0 0], 'EdgeColor', [0.8 0.8 0.8]);
                set(h, 'FaceAlpha', 0.1)
            end
            
            % === PLOT heat map of song dat - i.e. make sure they are
            % aligned - first extract min song dur here
            minsongsamps=min([minsongsamps length(segextract_tmp(j).songdat)]);
        end
        axis tight
        % ---- plot heat map of song amplitudes
        figure(hfig2);
        AllSongAmpl=nan(numtrials, minsongsamps);
        for j=1:numtrials
            [smdat, tt] = lt_plot_AmplSm(segextract_tmp(j).songdat(1:minsongsamps), fs);
            AllSongAmpl(j, :) = smdat;
        end
        %         subplotinds=sort([25 26]);
        hsplot = lt_subplot(2, 1, 1); hold on; hsplots =[hsplots hsplot];
        imagesc(tt, 1:numtrials, (AllSongAmpl));
        % ---- SMOOTHED MEAN FR
        figure(hfig1);
        hsplot=lt_subplot(8,4,29:30); hold on; hsplots=[hsplots hsplot];
        Yspks={};
        for j=1:numtrials
            inds=segextract_tmp(j).spk_Clust==clustnum;
            spktimes=segextract_tmp(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        % -- convert to smoothed rate
        [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
        % --- plot
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
        lt_plot_zeroline
        ylim([0 300]);
        clear segextract_tmp
        
        
        
        % ================================================================================== WN NO HIT TRIALS
        % (BOTH ESCAPE AND CATCH, CURRENTLY)
        segextract_tmp=segextract([segextract.hit_WN]~=1);
        
        % ===== SONG SPECTROGRAM
        hsplot= lt_subplot(8, 4, 3:4); hold on; hsplots =[hsplots hsplot];
        % then pick a random song
        tmp=randi(length(segextract_tmp),1);
        songdat=segextract_tmp(tmp).songdat;
        fs=segextract_tmp(tmp).fs;
        % --- plot
        lt_plot_spectrogram(songdat, fs, 1, 0);
        %        line([motif_predur motif_predur], ylim, 'Color','w');
        axis tight
        
        % ==== RASTER
        subplotinds=sort([7:4:23 8:4:25]);
        hsplot = lt_subplot(8, 4, subplotinds); hold on; hsplots =[hsplots hsplot];
        hsplots2=[hsplots2 hsplot];
        
        %         yind=1; % will increment to plot all
        
        numtrials=length(segextract_tmp);
        minsongsamps=[];
        for j=1:numtrials
            inds=segextract_tmp(j).spk_Clust==clustnum;
            spktimes=segextract_tmp(j).(spktimefield)(inds);
            
            if ~isempty(spktimes)
                for jj=1:length(spktimes)
                    spk=spktimes(jj);
                    line([spk spk], -[j-0.4 j+0.4], 'Color', plotcols{i}, 'LineWidth', 1);
                end
            end
            
            % ==== plot time of trial
            if plotTimesOnRaster==1
                trialtime=segextract_tmp(j).song_datenum;
                trialtime=datestr(trialtime, 'HH:MM');
                lt_plot_text(0, -j, trialtime, [0.4 0.4 0.4]);
            end
            
            % ========= WAS TRIAL A HIT?
            if segextract_tmp(j).hit_WN==1
                %                 line([-0.02 -0.002], [-j -j], 'Color','k', 'LineWidth' , 2);
                % note when WN hit
                X=[segextract_tmp(j).WNonset_sec  segextract_tmp(j).WNoffset_sec ...
                    segextract_tmp(j).WNoffset_sec  segextract_tmp(j).WNonset_sec];
                Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
                h=patch(X, Y, 'w', 'FaceColor', [0 0 0], 'EdgeColor', [0.8 0.8 0.8]);
                set(h, 'FaceAlpha', 0.1)
                
                % === PLOT heat map of song dat - i.e. make sure they are
                % aligned - first extract min song dur here
            end
            minsongsamps=min([minsongsamps length(segextract_tmp(j).songdat)]);
            
        end
        axis tight
        % ---- plot heat map of song amplitudes
        figure(hfig2);
        AllSongAmpl=nan(numtrials, minsongsamps);
        for j=1:numtrials
            [smdat, tt] = lt_plot_AmplSm(segextract_tmp(j).songdat(1:minsongsamps), fs);
            AllSongAmpl(j, :) = smdat;
        end
        hsplot = lt_subplot(2, 1, 2); hold on; hsplots =[hsplots hsplot];
        imagesc(tt, 1:numtrials, AllSongAmpl);
        
        % ---- SMOOTHED MEAN FR
        figure(hfig1);
        hsplot=lt_subplot(8,4, 31:32); hold on; hsplots=[hsplots hsplot];
        numtrials=length(segextract_tmp);
        Yspks={};
        for j=1:numtrials
            inds=segextract_tmp(j).spk_Clust==clustnum;
            spktimes=segextract_tmp(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        % -- convert to smoothed rate
        [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
        % --- plot
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
        lt_plot_zeroline
        ylim([0 300]);
        % -- overlay on the WN on plot
        lt_subplot(8,4,29:30);
        shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', [0.7 0.7 0.7]}, 1);
        
        % -- line for WN onset if available
        if ~isempty(WNchangeDateStrings)
            dstring=WNchangeDateStrings{1};
            dnum=datenum(dstring, 'ddmmmyyyy-HHMM');
            wnonind=sum([segextract_tmp.song_datenum]<dnum);
            line(xlim, -[wnonind+0.5 wnonind+0.5], 'Color','k');
        end
        
        
        clear segextract_tmp
        % === all subplots
        linkaxes(hsplots, 'x');
        linkaxes(hsplots2, 'xy');
        lt_subtitle(['neuron ' num2str(i) ', motif: ' motif_regexpr_str{m}]);
        
        
        
        % ========================= BINNED COMPARISONS OF SMOOTH RATE
        lt_figure; hold on;
        hsplots=[];
        hsplots2=[];
        segextract_WNon=segextract([segextract.hit_WN]==1);
        segextract_WNoff=segextract([segextract.hit_WN]~=1);
        buffertime=10; % 10 minutes pre and post limits determined by earlies and latest WN on song (to find WN off songs)
        
        % --- takes bin on N trials with WN on, and compares to all WN off
        % trials that are within that window +/- T minutes.
        N=NumTrialsBin;
        numtrials=length(segextract_WNon);
        assert(numtrials/N<=16, 'not enough subplot slots!!');
        for j=1:ceil(numtrials/N)
            trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
            hsplot=lt_subplot(5, 3, j); hold on;
            hsplots=[hsplots hsplot];
            hsplots2=[hsplots2 hsplot];
            
            Yspks={};
            for jj=trialsToPlot
                inds=segextract_WNon(jj).spk_Clust==clustnum;
                spktimes=segextract_WNon(jj).(spktimefield)(inds);
                Yspks=[Yspks spktimes];
            end
            % -- convert to smoothed rate
            [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            ytmp_text=max(ymean_hz);
            lt_plot_text(0, ytmp_text+40, ['WN on trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'r');
            
%             % -- save first bin to overlay over other bins
%             if j==1;
%                 xbin1=xbin;
%                 ymean1=ymean_hz;
%                 ysem1=ysem_hz;
%             else
%                 shadedErrorBar(xbin1, ymean1, ysem1, {'Color', [0.75 0.75 0.75]}, 1);
%             end

            % --- plot
            if length(trialsToPlot)>1 % mean
                shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{i}}, 1);
                lt_plot_zeroline
                ylim([0 300]);
            end
            
            % --- get WN OFF TRIALS WITHIN THIS TIME PERIOD
            periodON=segextract_WNon(trialsToPlot(1)).song_datenum;
            periodOFF=segextract_WNon(trialsToPlot(end)).song_datenum;
            % add buffer pre and post
            buffertime_days=buffertime/(24*60);
            periodON=periodON-buffertime_days;
            periodOFF=periodOFF+buffertime_days;
            trialsToPlot=find([segextract_WNoff.song_datenum] >= periodON ...
                & [segextract_WNoff.song_datenum] <= periodOFF);
            
            lt_plot_text(0, ytmp_text+20, ['WN off trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'k');
            
            Yspks={};
            for jj=trialsToPlot
                inds=segextract_WNoff(jj).spk_Clust==clustnum;
                spktimes=segextract_WNoff(jj).(spktimefield)(inds);
                Yspks=[Yspks spktimes];
            end
            % -- convert to smoothed rate
            [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            
            % --- plot
            if length(trialsToPlot)>1 % mean
                shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', 'k'}, 1);
            end
            
            % ----------------
            
            
            axis tight
        end
        
        linkaxes(hsplots, 'xy');
        
        % ================================ ANOTHER FIGURE; PLOT FF
        % TRAJECTORY
        FFvals=[segextract.FF_val];
        Tvals=[segextract.song_datenum];
        if ~isempty(FFvals)
            lt_figure; hold on;
            lt_subplot(1,2,1); hold on; % by rendition
            xlabel('rendition');
            plot(FFvals, 'o', 'Color', plotcols{i});
            % WN change points
            for dd=1:length(WNchangeDateStrings)
                datestring=WNchangeDateStrings{dd};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                allsongdnums=[segextract.song_datenum];
                sylind=find(allsongdnums>dnum, 1, 'first');
                
                line([sylind-0.5 sylind-0.5], ylim, 'Color', 'k');
                lt_plot_text(sylind+0.5, FFvals(1), datestring, 'k');
            end
            
            lt_subplot(1,2,2); hold on; % by time
            xlabel('time');
            plot(Tvals([segextract.hit_WN]==1), FFvals([segextract.hit_WN]==1), 'o', 'Color', plotcols{i}, 'MarkerFaceColor', plotcols{i});
            plot(Tvals([segextract.hit_WN]==0), FFvals([segextract.hit_WN]==0), 'o', 'Color', plotcols{i});
            datetick('x', 'ddmmm-HHMM', 'keepticks');
            % WN change points
            for dd=1:length(WNchangeDateStrings)
                datestring=WNchangeDateStrings{dd};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                line([dnum dnum], ylim, 'Color', 'k');
                %                 lt_plot_text(sylind+0.5, FFvals(1), datestring, 'k');
            end
            lt_subtitle(motif_regexpr_str{m})
            
        end
        
    end
    %
    %     pause; close all;
    %
end
