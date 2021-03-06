function lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong)

% plotIndivRaster; % one raster for each neuron/motif
% plotCombRast; % one figure, all rasters
% plotSmFR; % all smoothed FR.


%
% BirdToPlot = 'pu69wh78';
% % % ---- give it either
% % A) one neuron and a bunch of motifs or
% % B) bunch of neurons and one motif
% NeurToPlot = 1; % 4 % vector (e.g. [5 7]) - if [] then plots all;
% motiflist = {'a(b)', 'jbh(h)g'};
% plotbytime = 0; % links rasters for all motifs by time of song.

%%
% motifpredur = 0.15;
% motifpostdur = 0.15;

%%

AllRasters = {};
AllSmoothFR = {};
AllSmoothFR_sem = {};
AllSmoothFR_x = {};
AllNeurNum = [];
AllMotif = {};

numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    if ~strcmp(birdname, BirdToPlot)
        continue
    end
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        
        if ~isempty(NeurToPlot)
            if all(ii~=NeurToPlot)
                continue
            end
        end
        
        % ==================== PLOT
        for j = 1:length(motiflist)
            motiftoplot = motiflist{j};
            
            if plotIndivRaster==1
                lt_figure; hold on;
                hsplots = [];
            end
            
            % ================
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
            
            if isempty(SongDat.AllLabels)
                continue
            end
            
            collectWNhit = 0;
            preAndPostDurRelSameTimept = 1;
            RemoveIfTooLongGapDur = 1;
            FFparams.collectFF=0;
            clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                motiftoplot, motifpredur, motifpostdur, 1, '', FFparams, ...
                0, 1, collectWNhit, 0, 0, preAndPostDurRelSameTimept, RemoveIfTooLongGapDur, ...
                clustnum, PlotDirSong);
            
            if isempty(SegmentsExtract)
                continue
            end
            
            % =============================== DIR OR UNDIR?
            if PlotDirSong==0
                % then remove any dir
                SegmentsExtract([SegmentsExtract.DirSong]==1) = [];
            elseif PlotDirSong==1
                % then remove any undir
                SegmentsExtract([SegmentsExtract.DirSong]==0) = [];
            end
            
                
            
            % ================================
            
            
            % ------------- 1) PLOT RASTER
            if plotIndivRaster==1
                hsplot = lt_subplot(6,1,2:4); hold on;
                hsplots = [hsplots hsplot];
                title([birdname '-n' num2str(ii) '-' motiftoplot]);
                ylabel('trial, down is later');
                
                
                if plotbytime==1
                    for tt = 1:length(SegmentsExtract)
                        spktimes = SegmentsExtract(tt).spk_Times;
                        %             spktimes = spktimes(spktimes > WindowToPlot2(1) & ...
                        %                 spktimes < WindowToPlot2(2));
                        rendtime = SegmentsExtract(tt).global_tokenind_DatAlignedToOnsetOfThis;
                        for ttt =1:length(spktimes)
                            
                            line([spktimes(ttt) spktimes(ttt)], -[rendtime-10 rendtime+10], ...
                                'Color', 'k', 'LineWidth', 3);
                        end
                    end
                    axis tight
                else
                    for tt = 1:length(SegmentsExtract)
                        spktimes = SegmentsExtract(tt).spk_Times;
                        
                        %             spktimes = spktimes(spktimes > WindowToPlot2(1) & ...
                        %                 spktimes < WindowToPlot2(2));
                        for ttt =1:length(spktimes)
                            
                            line([spktimes(ttt) spktimes(ttt)], -[tt-0.4 tt+0.4], ...
                                'Color', 'k', 'LineWidth', 1);
                        end
                    end
                    axis tight
                    line([motifpredur motifpredur], ylim);
                    set(gca, 'Ytick', []);
                end
            end
            
            % ############################### COLLECT SPIKE TIMES
            AllRasters = [AllRasters {{SegmentsExtract.spk_Times}}];
            
            
            % ------------- 2) PLOT SMOOTHED FR
            SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, '');
            FRmat = [SegmentsExtract.FRsmooth_rate_CommonTrialDur];
            X = SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
            % 1) each trial
            if plotIndivRaster==1
                hsplot = lt_subplot(6,1,5); hold on;
                hsplots = [hsplots hsplot];
                for tt =1:length(SegmentsExtract)
                    plot(X, FRmat(:,tt), 'Color', [0.7 0.7 0.7]);
                end
                line([motifpredur motifpredur], ylim);
            end
            
            % 2) mean
            if length(SegmentsExtract)>2
                FRmean = mean(FRmat,2);
                FRsem = lt_sem(FRmat');
                
                if plotIndivRaster==1
                    hsplot = lt_subplot(6,1,6); hold on;
                    hsplots = [hsplots hsplot];
                    shadedErrorBar(X, FRmean, FRsem, {'Color','r'},1);
                    line([motifpredur motifpredur], ylim);
                end
                
                % ############################### COLLECT SMOOTHED FR
                AllSmoothFR = [AllSmoothFR FRmean];
                AllSmoothFR_sem = [AllSmoothFR_sem FRsem];
                AllSmoothFR_x = [AllSmoothFR_x X];
            end
            
            
            % ------------- 3) PLOT syl onset/offsets
            if plotIndivRaster==1
                hsplot = lt_subplot(6,1,1); hold on;
                %         hsplots = [hsplots hsplot];
                numsyltrialstoplot = 15;
                maxdur = motifpostdur+motifpredur -0.005;
                [SylContours, x] = lt_neural_v2_ANALY_GetSylContours(SegmentsExtract, ...
                    numsyltrialstoplot, maxdur);
                spy(SylContours);
                %         plot(x, SylContours, '+');
                set(gca, 'XTick', 1:50:size(SylContours,2), 'XTickLabel', 100*x(1:50:end))
                line([motifpredur motifpredur], ylim);
                
                % ---
                linkaxes(hsplots, 'x');
            end
            
            % ###################### collect stuff
            AllMotif = [AllMotif motiftoplot];
            AllNeurNum = [AllNeurNum ii];
        end
    end
end

%%

assert(length(AllSmoothFR) == length(AllNeurNum), 'asdfas');

%% ================= COMBINATION PLOT (RASTERS)
if plotCombRast ==1
    lt_figure; hold on;
    numplots = length(AllRasters);
    plotcols = lt_make_plot_colors(numplots, 0, 0);
    yval = 1;
    for i=1:numplots
        
        % ======= plot raster
        rasters = AllRasters{i};
        for j=1:length(rasters)
            spktimes = rasters{j};
            lt_neural_PLOT_rasterline(spktimes, yval, plotcols{i});
            yval = yval+1;
        end
        lt_plot_text(max(spktimes), yval, ['n' num2str(AllNeurNum(i)) '-' AllMotif{i}], plotcols{i});
    end
    axis tight
    ylabel('trial');
    xlabel('time (sec)');
    line([motifpredur motifpredur], ylim);
end

%% ================== COMBINATION PLOT (SMOOTHED FR)
if plotSmFR==1
    lt_figure; hold on;
    numplots = length(AllRasters);
    plotcols = lt_make_plot_colors(numplots, 0, 0);
    for i=1:numplots
        
        % ======= plot raster
        y = AllSmoothFR{i};
        ysem = AllSmoothFR_sem{i};
        x = AllSmoothFR_x{i};
        shadedErrorBar(x, y, ysem, {'Color', plotcols{i}},1);
        lt_plot_text(x(1), y(1), ['n' num2str(AllNeurNum(i)) '-' AllMotif{i}], plotcols{i});
    end
    axis tight
    xlabel('time (sec)');
    line([motifpredur motifpredur], ylim);
end


