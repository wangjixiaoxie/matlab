%% ########## take entier motif (stereotyped) and plot things. also compares RA and LMAN

% == TO DO: 
% trial by trial xcov.
% also subtract shuffle corrector.


%% ========== 1) extract full motifs and time warp

lt_neural_v2_ExtractFullMotifs;



%% === 2) plot all motifs, with smoothed firing, running std, for all cells
% === NOTE: will assert that has been timewarped (all neurons simultaneously) ...
i=1;
mm = 1;
plotcv = 0; % if 1, then plots running cv. if 0 then plots fr and running mean.

bregionsToPlot = {'RA', 'LMAN'};
plotcols = lt_make_plot_colors(length(bregionsToPlot), 0,0);
motifthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{mm};
% ------ one figure for each motif
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);

% ------- to collect
AllBregion = {};
AllFRsm = {};
AllFRsmActual = {}; % will be used if cv is taking up AllFRsm
AllOnsets = [];
AllOffsets =[];
for ii=1:numneurons
    
    bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).NOTE_Location;
    pcol = plotcols{strcmp(bregionsToPlot, bregionthis)};
    
    % ============= PLOT ALL TRIALS FOR THIS NEURON
    segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
    % -- make sure has been timewarped
    assert(isfield(segextract, 'LinTWSeg_OriginalSpkTimes'), 'need to time warp!!');
    
    segextract = lt_neural_SmoothFR(segextract);
    
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['n' num2str(ii) '-' motifthis '-'  bregionthis]);
    
    if plotcv==0
    plot(x, frmat, '-', 'Color', [0.7 0.7 0.7]);
    y = mean(frmat,2);
    plot(x, y, '-', 'LineWidth', 3, 'Color', pcol);
    elseif plotcv==1;
        
        y = std(frmat, 0, 2);
        y = y./mean(frmat,2);
        plot(x, y, '-', 'LineWidth', 3, 'Color', pcol);
    end
    
    AllFRsmActual = [AllFRsmActual; mean(frmat,2)];
        
    % --- overlay syl onsets/offsets
    onsets = segextract(1).motifsylOnsets;
    offsets = segextract(1).motifsylOffsets;
    
    axis tight
    lt_plot_zeroline;
    YLIM = ylim;
    
    for j=1:length(onsets);
        line([onsets(j) offsets(j)], [3 3], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-4 YLIM(2)-4], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'r');
        line([offsets(j) offsets(j)], ylim, 'Color', 'r');
    end
    
    
    % =========================== COLLECT FOR PLOTTING
    AllBregion = [AllBregion;  bregionthis];
    AllFRsm = [AllFRsm; y];
    AllOnsets = [AllOnsets; onsets];
    AllOffsets =[AllOffsets; offsets];
end

linkaxes(hsplots, 'x');


% ======================== PLOT OVELRAY
doztransform = 1; % fr trace (performs on the mean);
if plotcv==1
    doztransform=0;
end
lt_figure; hold on;

% ====================== 1) overlay all trials
lt_subplot(3,2,1); hold on;

tmp = cellfun('length', AllFRsm);
xmax = min(tmp);
for j=1:length(AllFRsm)
   
    pcol = plotcols{strcmp(bregionsToPlot, AllBregion{j})};
    
    y = AllFRsm{j}(1:xmax);
    x = 0:0.001:0.001*(length(y)-1);
    
    % -------- z-transform fr
    if doztransform==1
    y = zscore(y);
    end
    
    % ------- plot
    plot(x,y, '-', 'Color', pcol);
    
end

% ====================== 2) plot each with mean
splotcount = 1;
for i=1:length(bregionsToPlot)
    splotcount = splotcount+1;
    indstoplot = find(strcmp(AllBregion, bregionsToPlot{i}));
    lt_subplot(3,2,splotcount); hold on;
    title(bregionsToPlot{i});
    
    Yall = [];
    YallFR = [];
 for j=indstoplot'
   
    pcol = plotcols{strcmp(bregionsToPlot, AllBregion{j})};
    
    y = AllFRsm{j}(1:xmax);
    x = 0:0.001:0.001*(length(y)-1);
    
    % -------- z-transform fr
    if doztransform==1
    y = zscore(y);
    end
    
    % ------- plot
    plot(x,y, '-', 'Color', pcol);
    
    Yall = [Yall; y'];

    % ----- collect mean FR, even if doing cv
    YallFR = [YallFR; zscore(AllFRsmActual{j}(1:xmax))'];

 end

 % -- plot mean
 plot(x, mean(Yall,1), '-', 'LineWidth', 3, 'Color', pcol);
 
  
 if plotcv==1
     % -- then overlay fr mean
 plot(x, mean(YallFR,1), '-k');
 
 end

  % --- overlay syllables
    onsets = median(AllOnsets,1);
   offsets = median(AllOffsets,1);
   axis tight
    YLIM = ylim;

    for j=1:length(onsets);
        line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'r');
        line([offsets(j) offsets(j)], ylim, 'Color', 'r');
    end

end


% ====================== 2) overlay each bregion mean
    lt_subplot(3,2,splotcount+1); hold on;
    
for i=1:length(bregionsToPlot)
    
    indstoplot = find(strcmp(AllBregion, bregionsToPlot{i}));
    pcol = plotcols{i};
    Yall = [];
YallFR = [];
    for j=indstoplot'
   
    
    y = AllFRsm{j}(1:xmax);
    x = 0:0.001:0.001*(length(y)-1);
    
    % -------- z-transform fr
    if doztransform==1
    y = zscore(y);
    end
    
    Yall = [Yall; y'];
    
    % ----- collect mean FR, even if doing cv
    YallFR = [YallFR; zscore(AllFRsmActual{j}(1:xmax))'];
    
 end

 % -- plot mean
 shadedErrorBar(x, mean(Yall,1), lt_sem(Yall), {'Color', pcol},1);

end
 % --- overlay syllables
    onsets = median(AllOnsets,1);
   offsets = median(AllOffsets,1);
   axis tight
    YLIM = ylim;

    for j=1:length(onsets);
        line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'r');
        line([offsets(j) offsets(j)], ylim, 'Color', 'r');
    end
    
    

    
    % ====================== CROSS CORRELATION
    xcovwindmax = 0.1;
    if plotcv==0
        
        % --- only do if is for FR
        
        Ybyregion = [];
        % ============================== LMAN
        indstoplot = find(strcmp(AllBregion, 'LMAN'));
        Yall = [];
        for j=indstoplot'
            
            y = AllFRsm{j}(1:xmax);
            x = 0:0.001:0.001*(length(y)-1);
            
            % -------- z-transform fr
            if doztransform==1
                y = zscore(y);
            end
            
            Yall = [Yall; y'];
        end
        Ybyregion = [Ybyregion; mean(Yall,1)];
        
        
        % ============================== RA
        indstoplot = find(strcmp(AllBregion, 'RA'));
        Yall = [];
        for j=indstoplot'
            
            y = AllFRsm{j}(1:xmax);
            x = 0:0.001:0.001*(length(y)-1);
            
            % -------- z-transform fr
            if doztransform==1
                y = zscore(y);
            end
            
            Yall = [Yall; y'];
        end
        Ybyregion = [Ybyregion; mean(Yall,1)];
        
        
        % ================ DO XCORR
        [cc, lags] = xcov(Ybyregion(1,:), Ybyregion(2,:), xcovwindmax/0.001);
        lt_figure; hold on ;
        xlabel('LMAN ---- RA');
        ylabel('xcov');
        plot(lags*0.001, cc, '-k');
        lt_plot_zeroline_vert;
    end

