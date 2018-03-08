
%% ################################ DIR VS. UNDIR
% PLOT ALL MOTIFS - separating into DIR and UNDIR
% ORGANIZE MOTIFS BY SYLLABLE TYPE
close all; clear MOTIFSTATS_Compiled;
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 1;
saveOn = 0;
OrganizeByExpt =0;
collectFF=0;

% --- to make sure extracts motifs
% MotifsToCollect = {'pu69wh78', {'(j)jjbhhg', '(a)abhhg'}};
%     Params_regexp.motif_predur = 0.05;
%     Params_regexp.motif_postdur = 0.05;
%     Params_regexp.preAndPostDurRelSameTimept = 0;
%     Params_regexp.RemoveIfTooLongGapDur = 1;

MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF);


%% ====== PLOT ALL MOTIFS, OVERLAYING DIR AND UNDIR
close all;
i =1;
ii=3;
figcount=1;
subplotrows=4;
subplotcols=9;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str);
motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;
for j=1:nummotifs
    
    segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(j).SegmentsExtract;
    motifstr = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{j};
    
    % ==================== PLOT SMOOTHED FR (all trials + mean)
    % -- extract smooth fr
    clustnum = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).clustnum;
    segextract = lt_neural_SmoothFR(segextract, clustnum);
    
    % -------------- UNDIR
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([motifstr '-UNDIR']);
    inds = [segextract.DirSong]==0;
    
    frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
    y = mean(frmat,2);
    plot(x, y, '-k', 'LineWidth', 3);
    line([motifpredur motifpredur], ylim, 'Color','r');
    
    x1 = x;
    y1 = y;
    y1sem = lt_sem(frmat');
    
    
    % -------------- DIR
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('DIR');
    inds = [segextract.DirSong]==1;
    
    frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
    y = mean(frmat,2);
    plot(x, y, '-b', 'LineWidth', 3);
    line([motifpredur motifpredur], ylim, 'Color','r');
    
    x2 = x;
    y2 = y;
    y2sem = lt_sem(frmat');
    
    
    % -------------- COMBINED
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('UNDIR(bk) vs. DIR(bu)');
    
    shadedErrorBar(x1, y1, y1sem, {'Color', 'k'}, 1);
    shadedErrorBar(x2, y2, y2sem, {'Color', 'b'}, 1);
    line([motifpredur motifpredur], ylim, 'Color','r');
    
    
    % ---------------
    linkaxes(hsplots, 'xy');
    
end

%% =========== CONTEXT-DEPENDENCE REDUCED DURING DIR? [PLOT RAW]


i =1;
ii=3;
figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str;
motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;

% -- convert motifs to single syls.
[~, ~, ~, SingleSyls] = ...
    lt_neural_v2_extractSameType(motiflist, {motiflist{1}});
SingleSylsUnique = unique(SingleSyls);
clustnum = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).clustnum;

for j=1:length(SingleSylsUnique)
    singlesylthis = SingleSylsUnique{j};
    
    motifsthis = strfind(cell2mat(SingleSyls), singlesylthis);
    plotcols = lt_make_plot_colors(length(motifsthis), 0, 0);
    
    if length(motifsthis)<2
        continue
    end

    % ============================= UNDIR - OVERLAY ALL MOTIFS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['syl ' singlesylthis '-UNDIR']);
    
    for jj=1:length(motifsthis)
        mm = motifsthis(jj);
        
        segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
        motifstr = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{mm};
        segextract = lt_neural_SmoothFR(segextract, clustnum);
        
        % ================ extract smoothed fr for UNDIR and DIR
        inds = [segextract.DirSong]==0;
        
        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
        y = mean(frmat,2);
        ysem = lt_sem(frmat');
        
        shadedErrorBar(x, y, ysem, {'Color', plotcols{jj}},1);
    end
    %     axis tight
    line([motifpredur motifpredur], ylim, 'Color', 'k');
    
    
    % ============================= DIR - OVERLAY ALL MOTIFS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('DIR');
    
    for jj=1:length(motifsthis)
        mm = motifsthis(jj);
        
        segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
        motifstr = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{mm};
        segextract = lt_neural_SmoothFR(segextract, clustnum);
        
        % ================ extract smoothed fr for UNDIR and DIR
        inds = [segextract.DirSong]==1;
        
        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
        y = mean(frmat,2);
        ysem = lt_sem(frmat');
        
        shadedErrorBar(x, y, ysem, {'Color', plotcols{jj}},1);
    end
    %     axis tight
    line([motifpredur motifpredur], ylim, 'Color', 'k');
    
    
end

linkaxes(hsplots, 'xy');


%% =========== CONTEXT-DEPENDENCE REDUCED DURING DIR? [SUMMARY]

i =1;
ii=6;
windowpremotor = [-0.025 0.025]; % rel syl onset
equalizetrials = 1; % then subsamples to match trials between DIR and UNDIR

pairwisemetric = 'corr'; % corr of FR 1 vs 2
% pairwisemetric = 'absdiff'; % absolte diff in mean FR

% -------------------------------------
motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str;
motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;

% -- convert motifs to single syls.
[~, ~, ~, SingleSyls] = ...
    lt_neural_v2_extractSameType(motiflist, {motiflist{1}});
SingleSylsUnique = unique(SingleSyls);
clustnum = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).clustnum;

clear AllBranchDistances;
AllBranchDistances.UNDIR = [];
AllBranchDistances.DIR= [];

for j=1:length(SingleSylsUnique)
    singlesylthis = SingleSylsUnique{j};
    motifsthis = strfind(cell2mat(SingleSyls), singlesylthis);
    
    if length(motifsthis)<2
        continue
    end
    
    % ######################################################## UNDIR
    getdir =0;
    fname = 'UNDIR';
    
    % ----- COLLECT ALL FR WITHIN PREMOTOR WINDOW
    FRall = [];
    for jj=1:length(motifsthis)
        mm = motifsthis(jj);
        
        segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
        motifstr = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{mm};
        segextract = lt_neural_SmoothFR(segextract, clustnum);
        
        % ================ extract smoothed fr for UNDIR and DIR
        inds = find([segextract.DirSong]==getdir);
        
        if equalizetrials==1
        %- figure out min trials
        mintrials = min([sum([segextract.DirSong]==1) sum([segextract.DirSong]==0)]);
        % -- take mintrials number of trials from inds
        inds = randperm(length(inds), mintrials);        
        end
        
        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
        xinds = x>=windowpremotor(1)+motifpredur & x<windowpremotor(2)+motifpredur;
        
        frmat = frmat(xinds, :);
        
        if (1)
            % version 1, difference of means
        y = mean(frmat,2);
        FRall = [FRall; y'];
        
        else
            % version 2, mean of differences between all trials
%             ntrials = size(frmat,2);
        end
    end
      
    % ------- COMPUTE PAIRWISE DISTANCES 
    Yall =[];
    for jj=1:length(motifsthis)
        for jjj=jj+1:length(motifsthis)
        
            ythis = nan;
            if strcmp(pairwisemetric, 'absdiff')
            ythis = sum(abs(FRall(jj,:) - FRall(jjj,:)));
            elseif strcmp(pairwisemetric, 'corr');
                ythis = corr(FRall(jj,:)', FRall(jjj,:)');
            end
            
             Yall =[Yall ythis];           
        end
    end    
    AllBranchDistances.(fname) = [AllBranchDistances.(fname); mean(Yall)];
    
    
    
    % ######################################################## DIR
    getdir =1;
    fname = 'DIR';
    
    % ----- COLLECT ALL FR WITHIN PREMOTOR WINDOW
    FRall = [];
    for jj=1:length(motifsthis)
        mm = motifsthis(jj);
        
        segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
        motifstr = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{mm};
        segextract = lt_neural_SmoothFR(segextract, clustnum);
        
        % ================ extract smoothed fr for UNDIR and DIR
        inds = find([segextract.DirSong]==getdir);
        
        if equalizetrials==1
        %- figure out min trials
        mintrials = min([sum([segextract.DirSong]==1) sum([segextract.DirSong]==0)]);
        % -- take mintrials number of trials from inds
        inds = randperm(length(inds), mintrials);        
        end
        
        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
        xinds = x>=windowpremotor(1)+motifpredur & x<windowpremotor(2)+motifpredur;
        
        frmat = frmat(xinds, :);
        
        if (1)
            % version 1, difference of means
        y = mean(frmat,2);
        FRall = [FRall; y'];
        
        else
            % version 2, mean of differences between all trials
%             ntrials = size(frmat,2);
        end
    end
      
    % ------- COMPUTE PAIRWISE DISTANCES 
    Yall =[];
    for jj=1:length(motifsthis)
        for jjj=jj+1:length(motifsthis)
        
            ythis = nan;
            if strcmp(pairwisemetric, 'absdiff')
            ythis = sum(abs(FRall(jj,:) - FRall(jjj,:)));
            elseif strcmp(pairwisemetric, 'corr');
                ythis = corr(FRall(jj,:)', FRall(jjj,:)');
            end
            
             Yall =[Yall ythis];           
        end
    end    
    AllBranchDistances.(fname) = [AllBranchDistances.(fname); mean(Yall)];
    
end


% =============================== PLOT ALL PAIRS (I.E. ALL BRANCHES)
numbranch = length(AllBranchDistances.UNDIR);
lt_figure; hold on;
xlabel('UNDIR -- DIR');
ylabel('abs diff in FR across contexts');

x = [1 2];
y = [AllBranchDistances.UNDIR AllBranchDistances.DIR];
plot(x,y, '-ok');
xlim([0 3]);

for i=1:size(y,1)
    lt_plot_text(2.1, y(i, 2), ['motif' num2str(i)], 'r');
end