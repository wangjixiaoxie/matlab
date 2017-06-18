function lt_neural_v2_ContextFR(MOTIFSTATS_Compiled)


% kernelSD = []; % smoothing of FR
WindowToPlot = [-0.125 0.075]; % relative to onset of syl
numsyltrialstoplot = 5; % to overlay syl profiles, num to take for each neuron/motif.

%% ==== extract data for all syls across all motifs

% MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct);

%% ===== for each bird and expt, plot all contexts for each syl (separately for each neuron)
% ==== v1 - one plot for each motif (overlaying neurons)

numbirds = length(MOTIFSTATS_Compiled.birds);

for i=1:numbirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    for ii=1:numexpts
        tic
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        % ========= ONE SUBPLOT FOR EACH SYL/CONTEXT COMBINATION
        numsyls = length(MotifStats.params.motif_regexpr_str);
        numneurons = length(MotifStats.neurons);
        plotcols = lt_make_plot_colors(numneurons, 0, 0);
        
        figcount=1;
        subplotrows=7;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        
        % ---- go by order of single syls
        SingleSyls = MotifStats.params.singlesyls;
        plotcols_single = lt_make_plot_colors(length(SingleSyls), 0, 0);
        hsplots = [];
        ymax = [];
        for ss = 1:length(SingleSyls)
            singlesyl = SingleSyls{ss};
            
            % ---
            for j=1:numsyls
                
                sylname = MotifStats.params.motif_regexpr_str{j};
                
                % -- only continue if syl corresponds to this single syl
                syltmp = sylname(strfind(sylname, '(')+1);
                if ~strcmp(singlesyl, syltmp)
                    continue
                end
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([upper(syltmp) '-' sylname], 'Color', plotcols_single{ss});
                
                for nn=1:numneurons
                SylContours = []; % trial x timebin
                    
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    clustnum = MotifStats.neurons(nn).clustnum;
                    chan = MotifStats.neurons(nn).motif(j).Params.channel_board;
                    
                    if ~isfield(segextract, 'spk_Times')
                        continue
                    end
                    
                    % --- extract smoothed FR
                    segextract = lt_neural_SmoothFR(segextract, clustnum);
                    
                    % --- PLOT SMOOTHED FR
                    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    motif_predur = MotifStats.params.motif_predur;
                    inds = xbin>(motif_predur+WindowToPlot(1)) & xbin<(motif_predur+WindowToPlot(2));
                    
                    FRall = [segextract.FRsmooth_rate_CommonTrialDur];
                    
                    FRmean = mean(FRall(inds,:),2);
                    FRsem = lt_sem(FRall(inds,:)');
                    FRsem = FRsem';
                    X = xbin(inds);
                    
                    shadedErrorBar(X, FRmean, FRsem, {'Color', plotcols{nn}}, 1)
                    
                    lt_plot_text(0.9*X(end), 1.01*FRmean(end), ['ch' num2str(chan)], plotcols{nn});
                    ymax = max([ymax max(FRmean)]);
                    
                    % ----------- EXTRACT SUBSET OF ACOUSTIC ONSETS AND
                    % OFFSETS TO GET SENSE FOR SYL TIMINGS
                    sylconttmp = lt_neural_v2_ANALY_GetSylContours(segextract, ...
                        numsyltrialstoplot, WindowToPlot(2)+motif_predur);
                    SylContours = [SylContours; sylconttmp];
                    
                    % ==== OVERLAY SYL CONTOUR FOR THIS NEURON
                    x = (1:size(SylContours,2))/1000;
                    plot(x, -30+30*SylContours, 'Color', plotcols{nn});
                    %                     plot(segextract(trialindtmp).sylOnTimes_RelDataOnset, 1, 'og');
                    %                     plot(segextract(trialindtmp).sylOffTimes_RelDataOnset, 1, 'or');
%                     if sylname == 'd(d)d'; 
%                         keyboard
%                     end
                end
                % --
                try
                    line([motif_predur motif_predur], ylim, 'Color', 'k');
                    xlim([motif_predur+WindowToPlot(1) motif_predur+WindowToPlot(2)]);
                catch err
                end
            end
        end
        
        linkaxes(hsplots, 'xy');
        ylim([-40 ymax+10]);
        lt_subtitle([birdname '-' exptname]);
        toc
        pause;
    end
end


%% ===== for each bird and expt, plot all contexts for each syl (separately for each neuron)
% ==== v2 - one plot for each singlesyl/neuron, overlaying all motifs

numbirds = length(MOTIFSTATS_Compiled.birds);

for i=1:numbirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    for ii=1:numexpts
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        % ========= ONE SUBPLOT FOR EACH SYL/CONTEXT COMBINATION
        numsyls = length(MotifStats.params.motif_regexpr_str);
        numneurons = length(MotifStats.neurons);
        plotcols = lt_make_plot_colors(numneurons, 0, 0);
        
        figcount=1;
        subplotrows=7;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        
        % ---- go by order of single syls
        SingleSyls = MotifStats.params.singlesyls;
        plotcols_single = lt_make_plot_colors(length(SingleSyls), 0, 0);
        plotcols_motif = lt_make_plot_colors(numsyls, 0, 0);
        hsplots = [];
        ymax = [];
        for ss = 1:length(SingleSyls)
            singlesyl = SingleSyls{ss};
            
            for nn=1:numneurons
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                chan = MotifStats.neurons(nn).motif(1).Params.channel_board;
                batchname = MotifStats.neurons(nn).motif(1).Params.batchf;
                title([upper(singlesyl) '-ch' num2str(chan)], 'Color', plotcols_single{ss});
                
                % ---
                for j=1:numsyls
                    
                    SylContours = [];
                    sylname = MotifStats.params.motif_regexpr_str{j};
                    
                    % -- only continue if syl corresponds to this single syl
                    syltmp = sylname(strfind(sylname, '(')+1);
                    if ~strcmp(singlesyl, syltmp)
                        continue
                    end
                    
                    
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    clustnum = MotifStats.neurons(nn).clustnum;
                    chan = MotifStats.neurons(nn).motif(j).Params.channel_board;
                    
                    if ~isfield(segextract, 'spk_Times')
                        continue
                    end
                    
                    % --- extract smoothed FR
                    segextract = lt_neural_SmoothFR(segextract, clustnum);
                    
                    % --- PLOT SMOOTHED FR
                    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    motif_predur = MotifStats.params.motif_predur;
                    inds = xbin>(motif_predur+WindowToPlot(1)) & xbin<(motif_predur+WindowToPlot(2));
                    
                    FRall = [segextract.FRsmooth_rate_CommonTrialDur];
                    
                    FRmean = mean(FRall(inds,:),2);
                    FRsem = lt_sem(FRall(inds,:)');
                    FRsem = FRsem';
                    X = xbin(inds);
                    
                    shadedErrorBar(X, FRmean, FRsem, {'Color', plotcols_motif{j}}, 1)
                    
                    lt_plot_text(0.9*X(end), 1.01*FRmean(end), sylname, plotcols_motif{j});
                    ymax = max([ymax max(FRmean)]);
                    
                    % ==== OVERLAY SYL CONTOUR FOR THIS NEURON/MOTIF
                                        sylconttmp = lt_neural_v2_ANALY_GetSylContours(segextract, ...
                        numsyltrialstoplot, WindowToPlot(2)+motif_predur);
                                        SylContours = [SylContours; sylconttmp];

                    x = (1:size(SylContours,2))/1000;
                    plot(x, -30+30*SylContours, 'Color', plotcols_motif{j});
                end
                
                % --
                try
                line([motif_predur motif_predur], ylim, 'Color', 'k');
                xlim([motif_predur+WindowToPlot(1) motif_predur+WindowToPlot(2)]);
                catch  err
                end
            end
        end
        
        linkaxes(hsplots, 'xy');
        ylim([-40 ymax+10]);
        lt_subtitle([birdname '-' exptname '-' batchname]);
        
        pause;
    end
end





%%


