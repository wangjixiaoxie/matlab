function lt_neural_v2_ANALY_Swtch_Tcourse2(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg)

%%
% plotneurzscore = 1; then zscore; otherwise plots difference from base
plotFFbyRend=0;
%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);

WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

FFsmthbinsize = 10;

minrends = 5; % for both train and base

premotorWind = SwitchStruct.params.premotorWind;


% ------ for bins across renditions
minbinsize = 15; % will take this or numbase rends (divbided by 2), whichever is larger
maxbinsize = 25;
%%

% onlyPlotTargNontarg = 1; % if 1, then only targ syl, if 0, then all syls


%%


if mod(FFsmthbinsize,2)==0
    FFsmthbinsize= FFsmthbinsize+1; % conver to odd, so that median tval is at an actual datapoint.
end

%%
for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    if ~isempty(birdname_get)
        if ~strcmp(birdname, birdname_get)
            continue
        end
    end
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        if ~isempty(exptname_get)
            if ~strcmp(exptname, exptname_get)
                continue
            end
        end
        
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        targsyls = MotifStats.params.TargSyls;
        nummotifs = length(motiflist);
        
        WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
            MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            if ~isempty(switchnum_get)
                if switchnum_get ~= iii;
                    continue
                end
            end
            
            goodneurons = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons;
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
            swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            if onlyPlotTargNontarg==1
                motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls)) ...
                    find(ismember(motiflist, MotifStats.params.SameTypeSyls))];
            else
                motifstoplot = 1:nummotifs;
            end
            
            %             numcols = length(goodneurons)+1;
            
            % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
            for j=motifstoplot
                
                for nn=goodneurons
                    
                    lt_figure; hold on;
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    if ~isfield(segextract, 'fs')
                        continue
                    end
                    
                    %                     if length(baseInds)<minrends | length(trainInds) < minrends
                    %                         % even for good neurosn, could occur if some motifs
                    %                         % labeled pre but not post.
                    %                         continue
                    %                     end
                    
                    % =========================== FIGURE SOME THINGS OUT
                    % ABOUT BINSIZES
                    numrows = 3;
                    baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    %                     trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds_WithinDayOfSw);
                    trialstoplot = [baseInds trainInds];
                    
                    binsize = max([minbinsize length(baseInds)/2]);
                    binsize = min([binsize maxbinsize]);
                    
                    binOnsets = 1:binsize:length(trialstoplot);
                    binOffsets = [binOnsets(2:end)-1 length(trialstoplot)];
                    
                    binTrialInds = cell(1,length(binOnsets));
                    for kk=1:length(binOnsets)
                        binTrialInds{kk} = trialstoplot(binOnsets(kk):binOffsets(kk));
                    end
                    
                    
                    
                    % ######################################## A) plot FF
                    lt_subplot(1,numrows,1); hold on;
                    tvals = [segextract(trialstoplot).song_datenum];
                    ffvals = [segextract(trialstoplot).FF_val];
                    %
                    %
                    %                     tvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).tvals;
                    %                     ffvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).ffvals;
                    %
                    
                    % -- conver tvals to days
                    firstday = datestr(min(tvals), 'ddmmmyyyy');
                    tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
                    tvals = tvals.FinalValue;
                    swtimetmp = lt_convert_EventTimes_to_RelTimes(firstday, swthis);
                    swtimetmp = swtimetmp.FinalValue;
                    
                    % -- get zscore pitch if desired
                    if FFzscore==1
                    ffmean = mean(ffvals(tvals<swtimetmp));
                    ffstd = std(ffvals(tvals<swtimetmp));
                        
                    ffvals = (ffvals-ffmean)./ffstd;
                    end
                    
                    if plotFFbyRend==1
                        plot(ffvals, 1:length(ffvals), 'xk');
                    else
                        plot(ffvals, tvals, 'xk');
                    end
                    
                    
                    % --- plot means for each bin
                    tall = [];
                    ffall = [];
                    for kk =1:length(binTrialInds)
                        tmedian = median(tvals(binTrialInds{kk}));
                        tmean = mean(tvals(binTrialInds{kk}));
                        ffmean = mean(ffvals(binTrialInds{kk}));
                        ffsem = lt_sem(ffvals(binTrialInds{kk}));
                        
                        tall = [tall tmean];
                        ffall = [ffall ffmean];
                        
                        line([ffmean-ffsem ffmean+ffsem], [tmean tmean], 'Color','r', ...
                            'LineWidth', 2);
                    end
                    lt_plot(ffall, tall, {'Color','r', 'LineStyle', '-'});
                    
                    % ------------ line for mean baseline pitch
                    ffmean = mean(ffvals(tvals<swtimetmp));
                    ffstd = std(ffvals(tvals<swtimetmp));
                    line([ffmean ffmean], ylim, 'Color','k');
                    line([ffmean-ffstd ffmean-ffstd], ylim, 'Color', 'k', 'LineStyle', '--');
                    line([ffmean+ffstd ffmean+ffstd], ylim, 'Color', 'k', 'LineStyle', '--');
                    
                    % ----- line for baseline offset
                    line(xlim, [swtimetmp swtimetmp], 'Color','r');
                    
                    ylim([min(tvals)-0.02 max(tvals)+0.02])
                    ylabel('days');
                    
                    
                    
                    
                    
                    % ######################################## PLOT RASTER
                    % ---- extract spktimes
                    clustnum = MotifStats.neurons(nn).clustnum;
                    lt_subplot(1,numrows,2); hold on;
                    for tt = trialstoplot
                        indstmp = segextract(tt).spk_Clust==clustnum;
                        spktimes = segextract(tt).spk_Times(indstmp);
                        spktimes = spktimes(spktimes > WindowToPlot2(1) & ...
                            spktimes < WindowToPlot2(2));
                        lt_neural_PLOT_rasterline(spktimes, tt, 'k');
                    end
                    axis tight
                    set(gca, 'Ytick', []);
                    
                    % --- line for base vs. training
                    lastbaseind = max(baseInds);
                    line(xlim, -[-lastbaseind-0.5 -lastbaseind-0.5], 'Color','m', 'LineWidth', 2);
                    
                    % -- line for syl onset
                    line([MotifStats.params.motif_predur MotifStats.params.motif_predur], ...
                        ylim, 'Color', 'm', 'LineWidth', 2);
                    
                    % --- line for neural analysis window
                    line([MotifStats.params.motif_predur+premotorWind(1) ...
                        MotifStats.params.motif_predur+premotorWind(1)], ylim,  ...
                        'LineStyle' ,'-', 'Color', 'm');
                    
                    line([MotifStats.params.motif_predur+premotorWind(2) ...
                        MotifStats.params.motif_predur+premotorWind(2)], ylim,  ...
                        'LineStyle' ,'-', 'Color', 'm');
                    
                    % --- line dividing bins
                    for kk=1:length(binTrialInds)
                        triallast = binTrialInds{kk}(end);
                        line(xlim, -[-triallast-0.5 -triallast-0.5], 'Color','r');
                    end
                    
                    loct = SummaryStruct.birds(i).neurons(nn).NOTE_Location;
                    title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j} '-n' num2str(nn) '[' loct ']']);
                    
                    
                    % ######################################## PLOT SMOOTHED FR (BASE, TRAIN START,
                    % TRAIN END)
                    basefrtrials = [];
                    % ------------------ FOR EACH TIME BIN, PLOT
                    hsplots = [];
                    for kk = 1:length(binTrialInds)
                        inds = binTrialInds{kk};
                        
                        FRmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(inds(1)).FRsmooth_xbin_CommonTrialDur(:,1);
                        
                        frmean = mean(FRmat,2);
                        frsem = lt_sem(FRmat');
                        
                        % ---- collect if is baseline
                        if max(inds)<=max(baseInds)
                            % then this is a baseline bin
                            basefrtrials = [basefrtrials FRmat];
                        end
                        
                        % --- plot
                        plotpos = 3*length(binTrialInds)-3*(kk-1);
                        hsplot = lt_subplot(length(binTrialInds), 3, plotpos); hold on;
                        hsplots = [hsplots hsplot];
                        shadedErrorBar(x, frmean, frsem, {'Color','k'},1);
                        
                        % ------------ Separate by high and low FF
                        ffvals = [segextract(inds).FF_val];
                        if ~isempty(ffvals)
                            
                            ffmedian = median(ffvals);
                            
                            % --- high FF
                            frmat = FRmat(:, ffvals>=ffmedian);
                            frmean = mean(frmat,2);
                            frsem = lt_sem(frmat');
                            
                            shadedErrorBar(x, frmean, frsem, {'Color', 'r'},1);
                            
                            % --- low FF
                            frmat = FRmat(:, ffvals<ffmedian);
                            frmean = mean(frmat,2);
                            frsem = lt_sem(frmat');
                            
                            shadedErrorBar(x, frmean, frsem, {'Color', 'b'},1);
                        end
                        
                        % ---- overlay baseline, if this is not baseline
                        if max(inds)>max(baseInds)
                            % then at least part of this is after baseline
                            plot(x, mean(basefrtrials,2), '-w');
                        end
                        
                        
                        % -- line for syl onset
                        line([MotifStats.params.motif_predur MotifStats.params.motif_predur], ...
                            ylim, 'Color', 'm', 'LineWidth', 2);
                        
                        % --- line for neural analysis window
                        line([MotifStats.params.motif_predur+premotorWind(1) ...
                            MotifStats.params.motif_predur+premotorWind(1)], ylim,  ...
                            'LineStyle' ,'-', 'Color', 'm');
                        
                        line([MotifStats.params.motif_predur+premotorWind(2) ...
                            MotifStats.params.motif_predur+premotorWind(2)], ylim,  ...
                            'LineStyle' ,'-', 'Color', 'm');
                        
                    end
                    linkaxes(hsplots, 'xy');
                    YLIM = get(gca,'YLim');
                    ylim([0 YLIM(2)]);
                    xlim([0 MotifStats.params.motif_predur+MotifStats.params.motif_postdur]);
                    
                    
                    % =================
                    
                end
            end
        end
    end
end


