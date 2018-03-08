function lt_neural_POPLEARN_Plot(MOTIFSTATS_pop, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot, BregionWantedList, onlyPlotIfBothPrePostTrials)

DoUndir=1; % if 1, then only UNDIR. if 0, then only DIRected
% BregionWantedList = {{'LMAN', 'LMAN'}, {'LMAN', 'RA'}, {'RA', 'RA'}};
% onlyPlotIfBothPrePostTrials = 0; % if 1, then only plots neuronsets/switches
% for which neuron has data both pre and post switch.

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
                
                if isempty(inds_pre) & isempty(inds_post)
                    continue
                end
                
                if onlyPlotIfBothPrePostTrials ==1
                    if isempty(inds_pre) | isempty(inds_post)
                        continue
                    end
                end
                
                    
                        disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                
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
                    
                    % ---- separate dir and undir
                    isdir = [DAT.motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.DirSong]==1;
                    
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
                
                
                % ================= 2) CROSS CORR BASELINE, TRAIN FOR ALL SYLS
                figcount=1;
                subplotrows=4;
                subplotcols=8;
                fignums_alreadyused=[];
                hfigs=[];
                hsplots = [];
                
                for bbbb = 1:length(BregionWantedList)
%                     lt_figure; hold on;
                    bregionwanted = BregionWantedList{bbbb};
                    %                 bregionwanted = {'LMAN', 'RA'}; % will make sure xcov reflects this inputed order
                    
                    
                    for m=1:length(motiflist)
                        
                        motif = motiflist{m};
                        if isempty(DAT.motif(m))
                            continue
                        end
                        
                        if isempty(DAT.motif(m).XCov_neurpair)
                            continue
                        end
                        
                        % ====== find the neuron pairs that match desired brain
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
                        
                        
                        %% UNDIR
                        DoDIR = 0;
                        [figcount, subplotrows, subplotcols, fignums_alreadyused, hfigs, hsplots] = fn_plotxcov(DAT, m, indspairs, flippair, swthis, figcount, subplotrows, ...
                            subplotcols, fignums_alreadyused, hfigs, hsplots, motif, normToPSTH, DoDIR, bregionwanted);
                         %% DIR
                        DoDIR = 1;
                        [figcount, subplotrows, subplotcols, fignums_alreadyused, hfigs, hsplots] = fn_plotxcov(DAT, m, indspairs, flippair, swthis, figcount, subplotrows, ...
                            subplotcols, fignums_alreadyused, hfigs, hsplots, motif, normToPSTH, DoDIR, bregionwanted);
                       
                    end
                end
            end
            
        end
        
    end
end
end

function [figcount, subplotrows, subplotcols, fignums_alreadyused, hfigs, hsplots] = fn_plotxcov(DAT, m, indspairs, flippair, swthis, figcount, subplotrows, ...
    subplotcols, fignums_alreadyused, hfigs, hsplots, motif, normToPSTH, DoDIR, bregionwanted)

if DoDIR==1
    pcol = 'm';
else
    pcol = 'b';
end
if DoDIR==1
    if ~isfield(DAT.motif(m).XCov_neurpair(indspairs), 'DIR_ccRealAll')
        return
    end
    % -- then dirsong
    % ---- extract for those pairs
    ccRealAll = {DAT.motif(m).XCov_neurpair(indspairs).DIR_ccRealAll};
    ccShiftAll = {DAT.motif(m).XCov_neurpair(indspairs).DIR_ccShiftAll};
    ccPSTH = {DAT.motif(m).XCov_neurpair(indspairs).DIR_ccPSTH};
    xlags_sec = DAT.motif(m).XCov_neurpair(indspairs(1)).DIR_x;
else
    ccRealAll = {DAT.motif(m).XCov_neurpair(indspairs).ccRealAll};
    ccShiftAll = {DAT.motif(m).XCov_neurpair(indspairs).ccShiftAll};
    ccPSTH = {DAT.motif(m).XCov_neurpair(indspairs).ccPSTH};
    xlags_sec = DAT.motif(m).XCov_neurpair(indspairs(1)).x;
end


% -- flip any if needed
ccRealAll(flippair) = cellfun(@fliplr, ccRealAll(flippair), 'UniformOutput', 0);
ccShiftAll(flippair) = cellfun(@fliplr, ccShiftAll(flippair), 'UniformOutput', 0);
ccPSTH(flippair) = cellfun(@flipud, ccPSTH(flippair), 'UniformOutput', 0);


% ------------------ FIND BASE AND TRAINING INDS
indsundir = [DAT.motif(m).SegExtr_neurfakeID(1).SegmentsExtract.DirSong]==DoDIR;
songtimes = [DAT.motif(m).SegExtr_neurfakeID(1).SegmentsExtract(indsundir).song_datenum];

preInds = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
postInds = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
indtmp = length(postInds);
postInds_early = postInds(1:floor(indtmp/2));
postInds_late = postInds(floor(indtmp/2)+1:end);

assert(length(songtimes) == size(ccRealAll{1},1), 'asasdfasd');

%% PRE-SWITCH
if ~isempty(preInds)
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
        
        if length(preInds)>1
            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
        end
        plot(xlags_sec, ccnorm, 'Color', 'r');
        
        plot(xlags_sec, 10*ccfinal, '-', 'Color', pcol);
        
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    linkaxes(hsplots, 'xy');
end



%% ------------- POST (early)
if DoDIR==0
    % only split for UNDIR data (i.e. DIR does not have enough gnerally)
if ~isempty(postInds)
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
        
        if length(postInds_early)>1
            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
        end
        plot(xlags_sec, ccnorm, 'Color', 'r');
        
        plot(xlags_sec, 10*ccfinal, '-', 'Color', pcol);
        
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
        
        if length(postInds_late)>1
            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
        end
        plot(xlags_sec, ccnorm, 'Color', 'r');
        
        plot(xlags_sec, 10*ccfinal, '-', 'Color', pcol);
        
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    linkaxes(hsplots, 'xy');
end
else
        %% --------------- POST (late)
    % -- PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[POST(all)],' motif]);
    hsplots = [hsplots hsplot];
    if m==1
        ylabel('bk(raw), rd(shift), bu(minus,10x)');
    end
    
    for j=1:length(ccRealAll)
        
        % -- for each pair of neurons
        ccraw = mean(ccRealAll{j}(postInds,:), 1);
        ccraw_sem = lt_sem(ccRealAll{j}(postInds,:));
        ccshift = mean(ccShiftAll{j}(postInds,:), 1);
        ccshift_sem = lt_sem(ccShiftAll{j}(postInds,:));
        ccpst = ccPSTH{j}';
        if normToPSTH==1
            ccnorm = ccpst;
        else
            ccnorm = ccshift;
        end
        
        ccfinal = ccraw - ccnorm;
        
        if length(postInds)>1
            shadedErrorBar(xlags_sec, ccraw, ccraw_sem, {'Color', 'k'}, 1);
        end
        plot(xlags_sec, ccnorm, 'Color', 'r');
        
        plot(xlags_sec, 10*ccfinal, '-', 'Color', pcol);
        
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    linkaxes(hsplots, 'xy');
end


end
