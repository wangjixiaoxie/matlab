function lt_neural_v2_CTXT_FRanyclass(CLASSES, SummaryStruct, prms, plotPosControl, LMANorX)
%% lt 8/12/17 - plots fr for all classes of branch points

%%

numbirds = length(CLASSES.birds);
numsyltrialstoplot = 5;  % to overlay syl profiles, num to take for each neuron/motif.
WindowToPlot = []; % relative to onset of syl [ if empty, then plot all (and will minimize window to min dur]

Nmin = 5; % min dat to plot

%%



for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    birdname = CLASSES.birds(i).birdname;
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

    for ii=1:numneurons
        
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        if LMANorX==1
           if ~strcmp(SummaryStruct.birds(i).neurons(ii).NOTE_Location, 'LMAN')
               continue
           end
        elseif LMANorX==2
            if ~strcmp(SummaryStruct.birds(i).neurons(ii).NOTE_Location, 'X')
               continue
           end
        end
        
        for iii=1:numbranches
            
            % ========================== plot dat
            
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT')
                continue
            end
            
            if sum(cellfun('length', {CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum.SegmentsExtract})>5)<2
                % then fewer than 2 contexts with at least n =5;
                continue
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '-' CLASSES.birds(i).neurons(ii).branchnum(iii).regexprstr])

            numclasses = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum);
            plotcols = lt_make_plot_colors(numclasses, 0,0);
            ymax = [];
            xmin = [1000000000];
            
            for cc =1:numclasses
                
                
                segextract = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(cc).SegmentsExtract;
                sylname = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(cc).regexpstr;
                
                if ~isfield(segextract, 'fs')
                    continue
                end
                
                if numel(segextract)<Nmin
                    continue
                end
                
                % ---- extract shmoothed FR
                clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                segextract= lt_neural_SmoothFR(segextract, clustnum);
                
                % extract smoothed FR
                xtimes = segextract(1).FRsmooth_xbin_CommonTrialDur;
                if isempty(WindowToPlot)
                    inds = xtimes>0;
                else
                    inds = xtimes>=(prms.motifpredur+WindowToPlot(1)) & ...
                        xtimes<=(prms.motifpredur+WindowToPlot(2));
                end
                
                FRall = [segextract.FRsmooth_rate_CommonTrialDur];
                
                FRmean = mean(FRall(inds,:), 2);
                FRsem = lt_sem(FRall(inds,:)');
                xtimes = xtimes(inds);
                
                shadedErrorBar(xtimes, FRmean, FRsem, {'Color', plotcols{cc}},1);
                ymax = max([ymax max(FRmean)]);
                xmin = min([xmin xtimes(end)]);
                % -- annotate which motif
                lt_plot_text(xtimes(end), 1.01*FRmean(end), [sylname], plotcols{cc});
                
                % extract syl onsets/offsets
                sylconttmp = lt_neural_v2_ANALY_GetSylContours(segextract, ...
                    numsyltrialstoplot, xtimes(end-1));
                
                x = (1:size(sylconttmp,2))/1000;
                plot(x, -30+30*sylconttmp, 'Color', plotcols{cc});
            end
            
            axis tight;
            if ~isempty(ymax)
                ylim([-40 ymax+10]);
            end
            line([prms.motifpredur prms.motifpredur], ylim, 'Color', 'k');
            
            
            
            
            % ====================== PLOT POS CONTROL FOR THIS BRANCH, IF
            % EXISTS
            if plotPosControl==1
                
                
                if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT_POSCONTR')
                    continue
                end
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['[POSCON]' birdname '-n' num2str(ii) '-' CLASSES.birds(i).neurons(ii).branchnum(iii).regexprstr], 'Color', 'b')
                
                for cc =1:numclasses
                    
                    segextract = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract;
                    sylname = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(cc).regexpstr;
                    
                if ~isfield(segextract, 'fs')
                    continue
                end
                
                    
                    if numel(segextract)<Nmin
                        continue
                    end
                    
                                    % ---- extract shmoothed FR
                clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                segextract= lt_neural_SmoothFR(segextract, clustnum);

                    % extract smoothed FR
                    xtimes = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    if isempty(WindowToPlot)
                        inds = xtimes>0;
                    else
                        inds = xtimes>=(prms.motifpredur+WindowToPlot(1)) & ...
                            xtimes<=(prms.motifpredur+WindowToPlot(2));
                    end
                    
                    FRall = [segextract.FRsmooth_rate_CommonTrialDur];
                    
                    FRmean = mean(FRall(inds,:), 2);
                    FRsem = lt_sem(FRall(inds,:)');
                    xtimes = xtimes(inds);
                    
                    shadedErrorBar(xtimes, FRmean, FRsem, {'Color', plotcols{cc}},1);
                    ymax = max([ymax max(FRmean)]);
                    xmin = min([xmin xtimes(end)]);
                    % -- annotate which motif
                    lt_plot_text(xtimes(end), 1.01*FRmean(end), [sylname], plotcols{cc});
                    
                    % extract syl onsets/offsets
                    sylconttmp = lt_neural_v2_ANALY_GetSylContours(segextract, ...
                        numsyltrialstoplot, xtimes(end-1));
                    
                    x = (1:size(sylconttmp,2))/1000;
                    plot(x, -30+30*sylconttmp, 'Color', plotcols{cc});
                end
                
                axis tight;
                if ~isempty(ymax)
                    ylim([-40 ymax+10]);
                end
                line([prms.motifpredur prms.motifpredur], ylim, 'Color', 'k');
                
            end
            
        end
        
    end
    
end

