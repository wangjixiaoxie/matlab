function lt_neural_v2_DirUndir_LearnPlot(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot, plotIndTrial, SummaryStruct)
%% lt 2/25/18 - plots dir and undir smoothed FR, even if neuron doesn't
% have pre-switch data.

onlyPlotIfBothPrePostTrials = 0;
skipifnoDIR = 1; % skips motif

%% === plot all units

lt_neural_DISP_AllPopUnits(SummaryStruct);

%%
% BirdExptPairsToPlot = {};
% SwitchToPlot = [];
% plotIndTrial = 0;

%%
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ------------------------------- PLOT?
        if ~isempty(BirdExptPairsToPlot)
            
            ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
            ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
            
            if ~any(ind1+1 == ind2)
                disp(['SKIPPED ' birdname '-' exptname]);
                continue
            end
            
        end
        
        for iii=1:numswitches
            
            % ----------------------------------- PLOT?
            if ~isempty(SwitchToPlot)
                if ~any(SwitchToPlot == iii)
                    continue
                end
            end
            
            
            % ---- for this switch, figure out which neuron has data (at
            % laest one of pre or post
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
            numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
            motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
            motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;
            
            learningplotted = 0;
            for nn=1:numneurons
                bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).NOTE_Location;
                
                % ---------------- FIGURE OUT IF HAVE PRE AND POST SONG
                songtimes = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).Filedatenum_unsorted;
                
                inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                
                if isempty(inds_pre) & isempty(inds_post)
                    continue
                end
                
                if onlyPlotIfBothPrePostTrials ==1
                    if isempty(inds_pre) | isempty(inds_post)
                        continue
                    else
                        disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                    end
                end
                
                
                
                % ----------------------- PLOT EACH NEURON
                figcount=1;
                subplotrows=4;
                subplotcols=8;
                fignums_alreadyused=[];
                hfigs=[];
                hsplots =[];
                
                %                 ############################## PLOTS FOR THIS NEURON
                %                 =============== 1) PLOT LEARNING at target
                if learningplotted==0
                    lt_figure; hold on;
                    title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(nn)]);
                    targsyls = swthis.learningContingencies(1:2:end);
                    plotcolors = lt_make_plot_colors(length(targsyls), 0,0);
                    for m=1:length(targsyls)
                        tsyl = targsyls{m};
                        
                        indtmp = strcmp(motiflist, tsyl);
                        ffvals = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(indtmp).SegmentsExtract.FF_val];
                        tvals = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(indtmp).SegmentsExtract.song_datenum];
                        
                        % ---- separate dir and undir
                        isdir = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(indtmp).SegmentsExtract.DirSong]==1;
                        
                        % ----- UNDIR
                        plot(tvals(~isdir), ffvals(~isdir), 'o', 'Color', plotcolors{m});
                        % ---- DIR
                        plot(tvals(isdir), ffvals(isdir), 'o', 'Color', plotcolors{m}, 'MarkerFaceColor', 'k');
                        
                        learnconting = swthis.learningContingencies{2*m};
                        lt_plot_text(swthis.switchdnum, max(ffvals), [tsyl '-' num2str(learnconting)], plotcolors{m});
                    end
                    line([swthis.switchdnum swthis.switchdnum], ylim);
                    earliest = max([swthis.switchdnum_previous, floor(min(tvals))]);
                    latest = min([swthis.switchdnum_next, ceil(min(tvals))]);
                    
                    line([earliest earliest], ylim, 'LineStyle', '--');
                    line([latest latest], ylim, 'LineStyle', '--');
                    axis tight;
                    learningplotted = 1;
                end
                
                % ======================== 2) FOR EACH MOTIF, PLOT UNDIR
                % AND DIR
                for mm=1:length(motiflist)
                    
                    motifstr = motiflist{mm};
                    segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
                    
                    if isempty(segextract)
                        continue
                    end
                    
                    if skipifnoDIR==1
                        if ~any([segextract.DirSong])==1
                            continue
                        end
                    end
                    
                    % --- get pre and post inds
                    tvals = [segextract.song_datenum];
                    
                    indspre = tvals>swthis.switchdnum_previous ...
                        & tvals < swthis.switchdnum;
                    indspost = tvals>swthis.switchdnum ...
                        & tvals<swthis.switchdnum_next;
                    
                    % ---- get smoothed FR
                    segextract = lt_neural_SmoothFR(segextract);
                    
                    
                    %% ##################### PLOT PRE
                    indsEpoch = indspre;
                    EpochLab = 'PRE';
                    
                    if any(indsEpoch)
                        % -------------- UNDIR
                        inds = [segextract.DirSong]==0 & indsEpoch;
                        plottitle = [motifstr '-UNDIR'];
                        
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        
                        x1 = x;
                        y1 = y;
                        y1sem = lt_sem(frmat');
                        
                        
                        % -------------- DIR
                        inds = [segextract.DirSong]==1 & indsEpoch;
                        plottitle = [motifstr '-DIR'];
                        
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        x2 = x;
                        y2 = y;
                        y2sem = lt_sem(frmat');
                        
                        
                        % -------------- COMBINED
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title(motifstr);
                        xlabel(EpochLab);
                        
                        if mm==1
                            title([motifstr '-neur' num2str(nn)]);
                            ylabel([birdname '-' exptname '-sw' num2str(iii) '-' bregionthis]);
                        end
                        
                        shadedErrorBar(x1, y1, y1sem, {'Color', 'k'}, 1);
                        if size(frmat,2)>1
                            shadedErrorBar(x2, y2, y2sem, {'Color', 'b'}, 1);
                        end
                        line([motifpredur motifpredur], ylim, 'Color','r');
                    end
                    
                    %% ##################### PLOT POST
                    indsEpoch = indspost;
                    EpochLab = 'POST';
                    
                    if any(indsEpoch)
                        % -------------- UNDIR
                        inds = [segextract.DirSong]==0 & indsEpoch;
                        plottitle = [motifstr '-UNDIR'];
                        
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        
                        x1 = x;
                        y1 = y;
                        y1sem = lt_sem(frmat');
                        
                        
                        % -------------- DIR
                        inds = [segextract.DirSong]==1 & indsEpoch;
                        plottitle = [motifstr '-DIR'];
                        
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        
                        x2 = x;
                        y2 = y;
                        y2sem = lt_sem(frmat');
                        
                        
                        % -------------- COMBINED
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title(motifstr);
                        xlabel(EpochLab);
                        
                        if mm==1
                            title([motifstr '-neur' num2str(nn)]);
                            ylabel([birdname '-' exptname '-sw' num2str(iii) '-' bregionthis]);
                        end
                        
                        shadedErrorBar(x1, y1, y1sem, {'Color', 'k'}, 1);
                        if ~isempty(y2)
                            shadedErrorBar(x2, y2, y2sem, {'Color', 'b'}, 1);
                        end
                        line([motifpredur motifpredur], ylim, 'Color','r');
                    end
                end
                
                
                % ---------------
                if ~isempty(hsplots)
                    linkaxes(hsplots, 'xy');
                end
                
            end
            
        end
    end
end
