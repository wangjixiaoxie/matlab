function lt_neural_v2_ANALY_Learning(SummaryStruct, MOTIFSTATS, plottype, DivisorBaseSongs)
%% PARAMS

TargSyls = MOTIFSTATS.params.TargSyls;
motif_regexpr_str = MOTIFSTATS.params.motif_regexpr_str;
NumNeurons = length(MOTIFSTATS.neurons);
NumSyls = length(motif_regexpr_str);

motif_predur = MOTIFSTATS.params.motif_predur;
motif_postdur = MOTIFSTATS.params.motif_postdur;

%% Plot overview of experiment. (time, ffvals for targ and nontargs, and duration of recording for each neuron

% what are targ syls, same syl, diff syl

FirstDay = lt_neural_v2_ANALY_LearningPlot1(SummaryStruct, MOTIFSTATS);



%% FIGURE PARAMS
spktimefield = 'spk_Times';
PlotSmallTimeWindow = 1; % limit to window around syl onset
determineTrialBinBasedOnBaselineN =1;
% DivisorBaseSongs = 2;
BinsizeManual = 15; % num trials, if not using baseline N
MinNumRends = 10; % if fewer than this in baseline, then will take some from outside base.
window = 0.01; % neural;
windshift = 0.002;

%% FOR EACH NEURON, PLOT EACH MOTIF OVER TIME.
if ~strcmp(plottype, 'bysyl')
    for i=1:NumNeurons
        
        figcount=1;
        subplotrows=8;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        if strcmp(plottype, 'dotprod');
            subplotcols = 3; % add extra dot prod column.
        end
        
        WNchangeDateStrings = {};
        WNchangeDateStrings{1} = SummaryStruct.birds(1).neurons(i).LEARN_WNonDatestr;
        WNchangeDateStrings = [WNchangeDateStrings SummaryStruct.birds(1).neurons(i).LEARN_WNotherImportantDates];
        
        clustnum = MOTIFSTATS.neurons(i).clustnum;
        for ii=1:NumSyls
            
            sylname = motif_regexpr_str{ii};
            segextract = MOTIFSTATS.neurons(i).motif(ii).SegmentsExtract;
            numtrials = length(segextract);
            
            if ~isfield(segextract, 'song_datenum')
                continue
            end
            
            %  determine bins to use
            if determineTrialBinBasedOnBaselineN==1
                datestring=WNchangeDateStrings{1}; % assumes this is WN on.
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                N=ceil(sum([segextract.song_datenum] < dnum)/DivisorBaseSongs); % num, baseline rends divid by divisor.
                N = max([N MinNumRends]); % make sure at least this num of rends.
            else
                N=BinsizeManual;
            end
            
            
            % --- 1) FF
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            xlabel('rendition');
            ylabel('FF');
            title(['neuron ' num2str(i) ', ' sylname]);
            
            FFvals=[segextract.FF_val];
            Tvals=[segextract.song_datenum];
            if ~any(isnan(FFvals))
                
                %             lt_subplot(numrows, 3, ii*3-2); hold on; % by rendition
                
                tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
                Tvals_days = tmp.FinalValue;
                
                plot(Tvals_days, FFvals, 'o');
                axis tight;
                
                
                % WN change points IN PROGRESS
                %             for dd=1:length(WNchangeDateStrings)
                %                 datestring=WNchangeDateStrings{dd};
                %                 dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                %
                %                 allsongdnums=[segextract.song_datenum];
                %                 sylind=find(allsongdnums>dnum, 1, 'first');
                %
                %                 line(xlim, -[sylind-0.5 sylind-0.5], 'Color', 'k');
                %                 lt_plot_text(FFvals(1), -(sylind+0.5), datestring, 'k');
                %             end
                
                % ----- overlay +/- 1 SD of baseline
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                allsongdnums=[segextract.song_datenum];
                baseinds=allsongdnums<dnum;
                baseFF_mean = mean(FFvals(baseinds));
                baseFF_STD = std(FFvals(baseinds));
                
                line(xlim, [baseFF_mean baseFF_mean], 'Color', 'k');
                line(xlim, [baseFF_mean-baseFF_STD baseFF_mean-baseFF_STD], ...
                    'Color', 'k', 'LineStyle', '--');
                line(xlim, [baseFF_mean+baseFF_STD baseFF_mean+baseFF_STD], ...
                    'Color', 'k', 'LineStyle', '--');
                
                % ==== overlay mean FF for the windows used to extract FR
                RendMean = [];
                FFmeans = [];
                FFsems = [];
                FFstds = [];
                TimeMean = [];
                for j=1:ceil(numtrials/N);
                    trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                    
                    ffvalstmp = FFvals(trialsToPlot);
                    
                    ffmean = mean(ffvalstmp);
                    ffsem = lt_sem(ffvalstmp);
                    ffstd = std(ffvalstmp);
                    
                    TimeMean = [TimeMean mean(Tvals_days(trialsToPlot))];
                    RendMean = [RendMean mean(trialsToPlot)];
                    FFmeans = [FFmeans ffmean];
                    FFsems = [FFsems ffsem];
                    FFstds = [FFstds ffstd];
                end
                plot(TimeMean, FFmeans, 'o-', 'MarkerSize',8, 'Color', 'k', 'LineWidth', 2);
                for j=1:length(FFstds)
                    line([TimeMean(j) TimeMean(j)], [FFmeans(j)-FFstds(j) FFmeans(j)+FFstds(j)], ...
                        'Color', 'k', 'LineWidth', 2);
                end
                
                                if baseFF_mean-2*baseFF_STD < baseFF_mean+2*baseFF_STD
                ylim([baseFF_mean-2*baseFF_STD baseFF_mean+2*baseFF_STD]);
                end
 
                
                % -- line for base
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, dnum);
                WNontime = tmp.FinalValue;
                line([WNontime WNontime], ylim, 'Color','k');

            end
            
            
            % ---- 2) firing
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            plotcols = lt_make_plot_colors(ceil(numtrials/N), 1, [1 0 0]);
            for j=1:ceil(numtrials/N);
                trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                %             lt_plot_text(0, 250, ['trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'b');
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
                    if j==1
                        %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', 'b'}, 1);
                        plot(xbin, ymean_hz, 'Color', 'b', 'LineWidth', 2);
                    else
                        %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{j}}, 1);
                        plot(xbin, ymean_hz, '-', 'Color', plotcols{j}, 'LineWidth', 1.5);
                    end
                    lt_plot_zeroline
                    axis tight;
                    ylim([0 300]);
                    line([motif_predur motif_predur], ylim);
                end
            end
            
            % IN PROGRESS BELOW!!!!
            if (0)
                if strcmp(plottype, 'dotprod')
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    plotcols = lt_make_plot_colors(ceil(numtrials/N), 1, [1 0 0]);
                    for j=1:ceil(numtrials/N);
                        trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                        %             lt_plot_text(0, 250, ['trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'b');
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
                            if j==1
                                %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', 'b'}, 1);
                                plot(xbin, ymean_hz, 'Color', 'b', 'LineWidth', 2);
                            else
                                %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{j}}, 1);
                                plot(xbin, ymean_hz, '-', 'Color', plotcols{j}, 'LineWidth', 1.5);
                            end
                            lt_plot_zeroline
                            axis tight;
                            ylim([0 300]);
                            line([motif_predur motif_predur], ylim);
                        end
                    end
                    
                    
                end
            end
            
            
            
            % --- 3) spectrogram
            
        end
        
        lt_subtitle(['neuron ' num2str(i)]);
    end
end

%% ONE FIGURE FOR EACH MOTIF (ARRANGE ALL NEURONS)
if strcmp(plottype, 'bysyl')
    for ii=1:NumSyls
        sylname = motif_regexpr_str{ii};
        
        figcount=1;
        subplotrows=8;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        if strcmp(plottype, 'dotprod');
            subplotcols = 3; % add extra dot prod column.
        end
        
        neuroncount =1 ;
        plotcols_neurons = lt_make_plot_colors(NumNeurons, 0, 0);
        hsplots= [];
        for i=1:NumNeurons
            
            WNchangeDateStrings = {};
            WNchangeDateStrings{1} = SummaryStruct.birds(1).neurons(i).LEARN_WNonDatestr;
            WNchangeDateStrings = [WNchangeDateStrings SummaryStruct.birds(1).neurons(i).LEARN_WNotherImportantDates];
            
            clustnum = MOTIFSTATS.neurons(i).clustnum;
            
            segextract = MOTIFSTATS.neurons(i).motif(ii).SegmentsExtract;
            numtrials = length(segextract);
            
            if ~isfield(segextract, 'song_datenum')
                continue
            end
            
            %  determine bins to use
            if determineTrialBinBasedOnBaselineN==1
                datestring=WNchangeDateStrings{1}; % assumes this is WN on.
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                N=ceil(sum([segextract.song_datenum] < dnum)/DivisorBaseSongs); % num, baseline rends divid by divisor.
                N = max([N MinNumRends]); % make sure at least this num of rends.
            else
                N=BinsizeManual;
            end
            
            
            % --- 1) FF (OVERLAY FOR ALL NEURONS, SINCE SAME FF)
            if neuroncount==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            else
                lt_subplot(subplotrows, subplotcols, 1);
            end
            xlabel('rendition');
            ylabel('FF');
            title([sylname]);
            
            FFvals=[segextract.FF_val];
            Tvals=[segextract.song_datenum];
            if ~any(isnan(FFvals))
                
                %             lt_subplot(numrows, 3, ii*3-2); hold on; % by rendition
                
                tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
                Tvals_days = tmp.FinalValue;
                
                plot(Tvals_days, FFvals, 'o');
                axis tight;
                
                
                % WN change points IN PROGRESS
                %             for dd=1:length(WNchangeDateStrings)
                %                 datestring=WNchangeDateStrings{dd};
                %                 dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                %
                %                 allsongdnums=[segextract.song_datenum];
                %                 sylind=find(allsongdnums>dnum, 1, 'first');
                %
                %                 line(xlim, -[sylind-0.5 sylind-0.5], 'Color', 'k');
                %                 lt_plot_text(FFvals(1), -(sylind+0.5), datestring, 'k');
                %             end
                
                % ----- overlay +/- 1 SD of baseline
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                
                allsongdnums=[segextract.song_datenum];
                baseinds=allsongdnums<dnum;
                baseFF_mean = mean(FFvals(baseinds));
                baseFF_STD = std(FFvals(baseinds));
                
                line(xlim, [baseFF_mean baseFF_mean], 'Color', 'k');
                line(xlim, [baseFF_mean-baseFF_STD baseFF_mean-baseFF_STD], ...
                    'Color', 'k', 'LineStyle', '--');
                line(xlim, [baseFF_mean+baseFF_STD baseFF_mean+baseFF_STD], ...
                    'Color', 'k', 'LineStyle', '--');
                
                % ==== overlay mean FF for the windows used to extract FR
                RendMean = [];
                FFmeans = [];
                FFsems = [];
                FFstds = [];
                TimeMean = [];
                for j=1:ceil(numtrials/N);
                    trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                    
                    ffvalstmp = FFvals(trialsToPlot);
                    
                    ffmean = mean(ffvalstmp);
                    ffsem = lt_sem(ffvalstmp);
                    ffstd = std(ffvalstmp);
                    
                    TimeMean = [TimeMean mean(Tvals_days(trialsToPlot))];
                    RendMean = [RendMean mean(trialsToPlot)];
                    FFmeans = [FFmeans ffmean];
                    FFsems = [FFsems ffsem];
                    FFstds = [FFstds ffstd];
                end
                plot(TimeMean, FFmeans, 'o-', 'MarkerSize',8, 'Color', plotcols_neurons{i}, 'LineWidth', 2);
                for j=1:length(FFstds)
                    line([TimeMean(j) TimeMean(j)], [FFmeans(j)-FFstds(j) FFmeans(j)+FFstds(j)], ...
                        'Color', plotcols_neurons{i}, 'LineWidth', 2);
                end
                
                if baseFF_mean-2*baseFF_STD < baseFF_mean+2*baseFF_STD
                ylim([baseFF_mean-2*baseFF_STD baseFF_mean+2*baseFF_STD]);
                end
                
                % -- line for base
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                tmp = lt_convert_EventTimes_to_RelTimes(FirstDay, dnum);
                WNontime = tmp.FinalValue;
                line([WNontime WNontime], ylim, 'Color','k');
                
                % ----------- DRAW A LINE TELLING US WHERE WE GOT THIS NEURON
                % DATA FROM
                tmp = rand*0.5;
                line([min(Tvals_days) max(Tvals_days)], [baseFF_mean-(1.5+tmp)*baseFF_STD ...
                    baseFF_mean-(1.5+tmp)*baseFF_STD], 'Color', plotcols_neurons{i}, ...
                    'LineWidth', 1.5);
            end
            
            
            % ---- 2) firing
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            plotcols = lt_make_plot_colors(ceil(numtrials/N), 1, [1 0 0]);
            title(['neuron' num2str(i) '(chan' num2str(SummaryStruct.birds(1).neurons(i).channel) ')'], ...
                'Color', plotcols_neurons{i});
            for j=1:ceil(numtrials/N);
                trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                %             lt_plot_text(0, 250, ['trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'b');
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
                    if j==1
                        %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', 'b'}, 1);
                        plot(xbin, ymean_hz, 'Color', 'b', 'LineWidth', 2);
                    else
                        %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{j}}, 1);
                        plot(xbin, ymean_hz, '-', 'Color', plotcols{j}, 'LineWidth', 1.5);
                    end
                    lt_plot_zeroline
                    axis tight;
                    ylim([0 300]);
                    line([motif_predur motif_predur], ylim);
                end
            end
            
            % IN PROGRESS BELOW!!!!
            if (0)
                if strcmp(plottype, 'dotprod')
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    plotcols = lt_make_plot_colors(ceil(numtrials/N), 1, [1 0 0]);
                    for j=1:ceil(numtrials/N);
                        trialsToPlot=((j-1)*N+1):min(max(numtrials), j*N);
                        %             lt_plot_text(0, 250, ['trial ' num2str(trialsToPlot(1)) '-' num2str(trialsToPlot(end))], 'b');
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
                            if j==1
                                %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', 'b'}, 1);
                                plot(xbin, ymean_hz, 'Color', 'b', 'LineWidth', 2);
                            else
                                %                 shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{j}}, 1);
                                plot(xbin, ymean_hz, '-', 'Color', plotcols{j}, 'LineWidth', 1.5);
                            end
                            lt_plot_zeroline
                            axis tight;
                            ylim([0 300]);
                            line([motif_predur motif_predur], ylim);
                        end
                    end
                    
                    
                end
            end
            
            
            
            % --- 3) spectrogram
            neuroncount = neuroncount +1;
        end
        if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
        end
        %     lt_subtitle(['syl ' sylname]);
    end
end





