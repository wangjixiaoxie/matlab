function lt_neural_POPLEARN_Plot(MOTIFSTATS_pop, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot)
%% lt 118/18 - plots populaiton statistics relative to learning expts

% BirdExptPairsToPlot = {'pu69wh78', 'RALMANlearn1'}; % then only plots this expt

normToPSTH=0;
if normToPSTH==1
    disp('NOTE: psth not ready - subtracts cc of psth over entire recording day!');
    pause
end

%% PARAMS

%% RUN

numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        if ~isempty(BirdExptPairsToPlot)
           
            ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
            ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
            
            if ~any(ind1+1 == ind2)
                disp(['SKIPPED ' birdname '-' exptname]);
                continue
            end
            
        end
        
        for iii=1:numswitches
            
            if ~isempty(SwitchToPlot)
               if ~any(SwitchToPlot == iii)
                   continue
               end
            end
            
            % ---- for this switch, figure out which populations have data
            % overlapping the onset (i.e. has data both pre and post swictch)
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
            numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
            
            for ss = 1:numsets
                songfiles = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles{ss};
                songtimes = datenum(songfiles, 'yymmdd_HHMMSS');
                
                inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                
                if isempty(inds_pre) | isempty(inds_post)
                    continue
                else
                    disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                end
                
                
                
                % ############################################### ANALYSIS/PLOTS
                DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss);
                motiflist = {DAT.motif.regexpstr};
                
                % ================= 1) PLOT LEARNING AT TARGET(S)
                lt_figure; hold on;
                title([birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                targsyls = swthis.learningContingencies(1:2:end);
                plotcolors = lt_make_plot_colors(length(targsyls), 0,0);
                for m=1:length(targsyls)
                    tsyl = targsyls{m};
                    
                    indtmp = strcmp({DAT.motif.regexpstr}, tsyl);
                    tvals = [DAT.motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
                    ffvals = [DAT.motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.FF_val];
                    
                    plot(tvals, ffvals, 'o', 'Color', plotcolors{m});
                    learnconting = swthis.learningContingencies{2*m};
                    lt_plot_text(swthis.switchdnum, max(ffvals), [tsyl '-' num2str(learnconting)], plotcolors{m});
                end
                axis tight;
                line([swthis.switchdnum swthis.switchdnum], ylim);
                earliest = max([swthis.switchdnum_previous, floor(min(tvals))]);
                latest = min([swthis.switchdnum_next, ceil(min(tvals))]);
                
                line([earliest earliest], ylim, 'LineStyle', '--');
                line([latest latest], ylim, 'LineStyle', '--');
                
                
                % ================= 2) CROSS CORR BASELINE, TRAIN FOR ALL SYLS
                figcount=1;
                subplotrows=4;
                subplotcols=length(motiflist);
                fignums_alreadyused=[];
                hfigs=[];
                hsplots = [];
                
                BregionWantedList = {{'LMAN', 'LMAN'}, {'LMAN', 'RA'}, {'RA', 'RA'}};
                for bbbb = 1:length(BregionWantedList)
                lt_figure; hold on;
                    bregionwanted = BregionWantedList{bbbb};
                    %                 bregionwanted = {'LMAN', 'RA'}; % will make sure xcov reflects this inputed order
                    
                    
                    for m=1:length(motiflist)
                        motif = motiflist{m};
                        
                        % ---- find the neuron pairs that match desired brain
                        % regions - note if need to fliplr the xcov
                        bregionlist = {DAT.motif(m).XCov_neurpair(:).bregtmp};
                        indspairs = [];
                        flippair = [];
                        disp(['--------- DESIRED ' bregionwanted]);
                        for bb = 1:length(bregionlist)
                            if all(strcmp(bregionlist{bb}, bregionwanted))
                                % -- keep and don't flip
                                indspairs = [indspairs bb];
                                flippair = [flippair 0];
                            elseif all(strcmp(fliplr(bregionlist{bb}), bregionwanted))
                                % -- keep, and flip
                                indspairs = [indspairs bb];
                                flippair = [flippair 1];
                            end
                        end
                        flippair = logical(flippair);
                        
                        if isempty(indspairs)
                            % then no neuron pairs with these brain regions
                            continue
                        end
                        
                        % ---- extract for those pairs
                        ccRealAll = {DAT.motif(m).XCov_neurpair(indspairs).ccRealAll};
                        ccShiftAll = {DAT.motif(m).XCov_neurpair(indspairs).ccShiftAll};
                        ccPSTH = {DAT.motif(m).XCov_neurpair(indspairs).ccPSTH};
                        xlags_sec = DAT.motif(m).XCov_neurpair(indspairs(1)).x;
                        
                        % -- flip any if needed
                        ccRealAll(flippair) = cellfun(@fliplr, ccRealAll(flippair), 'UniformOutput', 0);
                        ccShiftAll(flippair) = cellfun(@fliplr, ccShiftAll(flippair), 'UniformOutput', 0);
                        ccPSTH(flippair) = cellfun(@flipud, ccPSTH(flippair), 'UniformOutput', 0);
                        
                        % ------------------ FIND BASE AND TRAINING INDS
                        songtimes = [DAT.motif(m).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
                        preInds = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                        postInds = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                        indtmp = length(postInds);
                        postInds_early = postInds(1:floor(indtmp/2));
                        postInds_late = postInds(floor(indtmp/2)+1:end);
                        
                        assert(length(songtimes) == size(ccRealAll{1},1), 'asasdfasd');
                        %% PRE-SWITCH
                        
                        % -- PLOT
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['[PRE],' motif]);
                        hsplots = [hsplots hsplot];
                        if m==1
                            ylabel('bk(raw), rd(shift), bu(minus,10x)');
                            xlabel([bregionwanted{1} '-' bregionwanted{2}])
                        end
                        
                        for j=1:length(ccRealAll)
                            
                            % -- for each pair of neurons
                            ccraw = mean(ccRealAll{j}(preInds,:), 1);
                            ccraw_sem = lt_sem(ccRealAll{j}(preInds,:));
                            ccshift = mean(ccShiftAll{j}(preInds,:), 1);
                            ccshift_sem = lt_sem(ccShiftAll{j}(preInds,:));
                            ccpst = ccPSTH{j}';
                            if normToPSTH==1
                                ccnorm = ccpst;
                            else
                                ccnorm = ccshift;
                            end
                            
                            ccfinal = ccraw - ccnorm;
                            
                            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
                            plot(xlags_sec, ccnorm, 'Color', 'r');
                            
                            plot(xlags_sec, 10*ccfinal, '-b');
                            
                        end
                        
                        axis tight;
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                        
                        
                        linkaxes(hsplots, 'xy');
                        
                        
                        %% ------------- POST (early)
                        % -- PLOT
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['[POST(early)],' motif]);
                        hsplots = [hsplots hsplot];
                        if m==1
                            ylabel('bk(raw), rd(shift), bu(minus,10x)');
                        end
                        
                        for j=1:length(ccRealAll)
                            
                            % -- for each pair of neurons
                            ccraw = mean(ccRealAll{j}(postInds_early,:), 1);
                            ccraw_sem = lt_sem(ccRealAll{j}(postInds_early,:));
                            ccshift = mean(ccShiftAll{j}(postInds_early,:), 1);
                            ccshift_sem = lt_sem(ccShiftAll{j}(postInds_early,:));
                            ccpst = ccPSTH{j}';
                            if normToPSTH==1
                                ccnorm = ccpst;
                            else
                                ccnorm = ccshift;
                            end
                            
                            ccfinal = ccraw - ccnorm;
                            
                            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
                            plot(xlags_sec, ccnorm, 'Color', 'r');
                            
                            plot(xlags_sec, 10*ccfinal, '-b');
                            
                        end
                        
                        axis tight;
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                        
                        
                        linkaxes(hsplots, 'xy');
                        
                        %% --------------- POST (late)
                        % -- PLOT
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['[POST(late)],' motif]);
                        hsplots = [hsplots hsplot];
                        if m==1
                            ylabel('bk(raw), rd(shift), bu(minus,10x)');
                        end
                        
                        for j=1:length(ccRealAll)
                            
                            % -- for each pair of neurons
                            ccraw = mean(ccRealAll{j}(postInds_late,:), 1);
                            ccraw_sem = lt_sem(ccRealAll{j}(postInds_late,:));
                            ccshift = mean(ccShiftAll{j}(postInds_late,:), 1);
                            ccshift_sem = lt_sem(ccShiftAll{j}(postInds_late,:));
                            ccpst = ccPSTH{j}';
                            if normToPSTH==1
                                ccnorm = ccpst;
                            else
                                ccnorm = ccshift;
                            end
                            
                            ccfinal = ccraw - ccnorm;
                            
                            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
                            plot(xlags_sec, ccnorm, 'Color', 'r');
                            
                            plot(xlags_sec, 10*ccfinal, '-b');
                            
                        end
                        
                        axis tight;
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                        
                        
                        linkaxes(hsplots, 'xy');
                        
                        
                        
                    end
                end
            end
            
        end
        
    end
end

