function lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime)




%%

numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    if ~strcmp(birdname, BirdToPlot)
        continue
    end
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        
        if ii~=NeurToPlot
            continue
        end
        
        
        % ==================== PLOT
        for j = 1:length(motiflist)
            motiftoplot = motiflist{j};
            
            lt_figure; hold on;
            hsplots = [];
            % ================
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
            
            
            motifpredur = 0.17;
            motifpostdur = 0.22;
            collectWNhit = 0;
            preAndPostDurRelSameTimept = 1;
            RemoveIfTooLongGapDur = 1;
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                motiftoplot, motifpredur, motifpostdur, 1, '', '', ...
                0, 1, collectWNhit, 0, 0, preAndPostDurRelSameTimept, RemoveIfTooLongGapDur);
            
            
            % ------------- 1) PLOT RASTER
            if plotbytime==1
            hsplot = lt_subplot(6,1,2:4); hold on;
            hsplots = [hsplots hsplot];
            title([birdname '-n' num2str(ii) '-' motiftoplot]);
            ylabel('trial, down is later');
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
            hsplot = lt_subplot(6,1,2:4); hold on;
            hsplots = [hsplots hsplot];
            title([birdname '-n' num2str(ii) '-' motiftoplot]);
            ylabel('trial, down is later');
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
            
            % ------------- 2) PLOT SMOOTHED FR
            SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, '');
            FRmat = [SegmentsExtract.FRsmooth_rate_CommonTrialDur];
            X = SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
            % 1) each trial
            hsplot = lt_subplot(6,1,5); hold on;
            hsplots = [hsplots hsplot];
            for tt =1:length(SegmentsExtract)
                plot(X, FRmat(:,tt), 'Color', [0.7 0.7 0.7]);
            end
                        line([motifpredur motifpredur], ylim);

                        % 2) mean
            if length(SegmentsExtract)>2
                hsplot = lt_subplot(6,1,6); hold on;
                hsplots = [hsplots hsplot];
                FRmean = mean(FRmat,2);
                FRsem = lt_sem(FRmat');
                shadedErrorBar(X, FRmean, FRsem, {'Color','r'},1);
            end
                        line([motifpredur motifpredur], ylim);

            
            % ------------- 3) PLOT syl onset/offsets
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
    end
end

