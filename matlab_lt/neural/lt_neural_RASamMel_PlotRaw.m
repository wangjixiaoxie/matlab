function lt_neural_RASamMel_PlotRaw(SummaryStruct)
%% Params

dirBirdsRaw = '/bluejay5/stork_RA_physiology/birds/';

%% Lt 9/24/17 - ==================== FOR EACH NEURON, PLOT
% 1) SONGS, LABELS; 2) RAW DAT; 3) RASTERS

%%

numbirds = length(SummaryStruct.birds);
count = 0;
lt_figure; hold on;
for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    numneurons = length(SummaryStruct.birds(i).neurons);
    for ii=1:numneurons
        
        datstruct = SummaryStruct.birds(i).neurons(ii).dat;
        filedate = SummaryStruct.birds(i).neurons(ii).date;
        
        pcterror = median(datstruct.pct_error);
        pctISIleq1ms = median(datstruct.pct_ISIs_leq_1ms);
        
        lt_figure; hold on;
        hsplots = [];
        % ================= 1) PLOT A RANDOM SONG FILE
        if (0) % should get directories and save hard copy to get with summary struct
            % NOTE: should also save dates of song filenames
            
            randsong = randi(length(datstruct.fname_arr));
            datstruct.fname_arr(randsong);
            
            cd(dirBirdsRaw)
            cd(birdname)
            
            %        tmp = eval(['!find . -name ' datstruct.fname_arr{randsong}])
            %        [~, tmp] = system(['find . -name ' datstruct.fname_arr{randsong}])
            %        [~, tmp] = system(['find . -name ' [datstruct.fname_arr{randsong} '.not.mat']]);
            [~, tmp] = system(['find . -name ' [datstruct.fname_arr{randsong} '.not.mat'] ' -print | head -n 1']);
            
            %        assert(size(tmp,1)==1, 'asdfasd');
            tmp = tmp(1,:);
            slashes = strfind(tmp, '/');
            cd(tmp(1:slashes(end)-1)); % file should be here
            disp(tmp(1:slashes(end)-1));
        end
        
        
        % =============== 2) extract a given syllable
        sylsunique = unique(datstruct.labels);
        randsyl = sylsunique(randi(length(sylsunique)));
        
        % -------- run regexp extract
        [SongDat, NeurDat, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, i, ii);

        motifpredur = 0.15;
        motifpostdur = 0.15;
        collectWNhit = 0;
        preAndPostDurRelSameTimept = 1;
        RemoveIfTooLongGapDur = 1;
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            ['(' randsyl ')'], motifpredur, motifpostdur, 1, '', '', ...
            0, 1, collectWNhit, 0, 0, preAndPostDurRelSameTimept, RemoveIfTooLongGapDur);
        
        
        % ------------- 1) PLOT RASTER
        hsplot = lt_subplot(6,1,2:4); hold on;
        hsplots = [hsplots hsplot];
        title(['pcterror:' num2str(pcterror) ', pctISI<=1ms:' num2str(pctISIleq1ms)]);
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
        set(gca, 'Ytick', []);
        
        
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
        % 2) mean
        if length(SegmentsExtract)>2
        hsplot = lt_subplot(6,1,6); hold on;
        hsplots = [hsplots hsplot];
        FRmean = mean(FRmat,2);
        FRsem = lt_sem(FRmat');
        shadedErrorBar(X, FRmean, FRsem, {'Color','r'},1);
        end
        
        % ------------- 3) PLOT syl onset/offsets
        hsplot = lt_subplot(6,1,1); hold on;
        title([birdname '-n' num2str(ii) '-' randsyl]);
%         hsplots = [hsplots hsplot];
        numsyltrialstoplot = 15;
        maxdur = motifpostdur+motifpredur -0.005;
        [SylContours, x] = lt_neural_v2_ANALY_GetSylContours(SegmentsExtract, ...
            numsyltrialstoplot, maxdur);
        spy(SylContours);
        %         plot(x, SylContours, '+');
        set(gca, 'XTick', 1:25:size(SylContours,2), 'XTickLabel', x(1:25:end))
        
        
        
        % ---
        linkaxes(hsplots, 'x');
        
        % ---- sanity check - units correct
        if (0)
            count = count+1;
        plot(NeurDat.spiketimes(end), count, 'ok');
        plot(NeurDat.TotalSamps./NeurDat.metaDat(1).fs, count, 'ok');
        plot(SongDat.AllOffsets(end), count, 'ok');
        tmptmp = [NeurDat.spiketimes(end) NeurDat.TotalSamps./NeurDat.metaDat(1).fs SongDat.AllOffsets(end)];
        line([min(tmptmp) max(tmptmp)], [count count], 'Color', 'r');
        title('all red lines should be small!! - means units are all matched');
        xlabel('rec duration (3 dots are spk, offsets, samps)');
        ylabel('neuron');
        end
        
        % -- 
        disp('PAUSED --- press keyboard');
        pause
        close all
    end
    
end