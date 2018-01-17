plotRaw = 0;
plotRaw_CC = 0;
onlySetsWithRALMANpaired=1;
tempfix=1;
handcodeinteraction = 1; % gets interaction variable by hand ..


plotCols = lt_make_plot_colors(3, 0, 0);
plotColsRegions = {'LMAN', 'RA', 'X'};


% ========== xcorr analysis of frate
xcovwindmax = 0.1;
xcov_dattotake = [-0.125 0.075]; % rel syl onset

%%
NumBirds = length(MOTIFSTATS_pop.birds);

% =============== SINGLE NEURON VS. FF
All_NeurFF_Rho = [];
All_NeurFF_P = [];
All_NeurFF_Brain = {};

% =============== NEURON PAIRS VS. FF
All_NeurPairFF_beta = []; % trials x 3(2 regions + interactio)
All_NeurPairFF_pval = [];
All_NeurPairFF_Brain = {};
All_NeurPairFF_NeurID = [];
All_NeurPairFF_BirdID = [];
All_NeurPairFF_ExptID = [];
All_NeurPairFF_Ntrials = [];
All_NeurPairFF_motifnum = [];
All_NeurPairFF_setnum = [];

% =============== BETWEEN NEURONS
All_NeurNeur_Rho = [];
All_NeurNeur_P = [];
All_NeurNeur_Brain = {};

% =============== BETWEEN NEURONS (XCORR)
All_NeurNeur_xcorr = {};
All_NeurNeur_xcorr_lags = [];
All_NeurNeur_xcorr_shift = {};
All_NeurNeur_FRreg1 = {};
All_NeurNeur_FRreg2 = {};

for i=1:NumBirds
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
    params = MOTIFSTATS_pop.birds(i).params;
    
    for ii=1:numexpts
        
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
        for iii=1:numsets
            disp(['bird' num2str(i) '-expt' num2str(ii) '-set' num2str(iii)])
            
            neurons_thisset = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{iii};
            nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif);
            %             DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii);
            nneur = length(neurons_thisset);
            
            % =========== WHAT BRAIN REGIONS ARE THESE?
            BrainregionsList = {SummaryStruct.birds(i).neurons(neurons_thisset).NOTE_Location};
            
            
            
            %% ######################### LMAN/RA PAIRS [MODIFY TO BE GENERAL]
            indsLMAN = strcmp(BrainregionsList, 'LMAN');
            indsRA = strcmp(BrainregionsList, 'RA');
            indsX = strcmp(BrainregionsList, 'X');
            
            if onlySetsWithRALMANpaired==1
                if ~any(indsLMAN) & ~any(indsRA)
                    continue
                end
            end
            
            %%
            for mm = 1:nummotifs
                
                segextract_for_trialdur = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                motifstr = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).regexpstr;
                motifpredur = MOTIFSTATS_pop.birds(i).params.motif_predur;
                
                
                %% ################## THINGS RELATED TO FF
                
                %                 if ~isfield(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm), 'POST')
                %                     continue
                %                 end
                %                 if isempty(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST)
                %                     % then no FF dat, skip
                %                     continue
                %                 end
                %
                %
                %                 % ============================= COLLECT DATA
                %                 FF = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.FF;
                %                 NspksAll = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.Nspks_allneur;
                
                
                if ~isfield(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm), 'POST')
                    FF = nan(1,nneur);
                    NspksAll = nan(1,nneur);
                elseif isempty(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST)
                    FF = nan(1,nneur);
                    NspksAll = nan(1,nneur);
                else
                    % ============================= COLLECT DATA
                    FF = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.FF;
                    NspksAll = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.Nspks_allneur;
                    
                end
                
                
                
                % ========== for plotting raw
                figcount=1;
                subplotrows=5;
                subplotcols=3;
                fignums_alreadyused=[];
                hfigs=[];
                
                
                %% 1) PERFORM MULTILINEAR REGRESSION [ALL PAIRS]
                for nn=1:nneur
                    for nnn=nn+1:nneur
                        
                        % --- brain regions
                        bregtmp = {BrainregionsList{nn}, BrainregionsList{nnn}};
                        
                        
                        
                        if isnan(FF(1))
                            % then skip this neuron pair
                            All_NeurPairFF_beta = [All_NeurPairFF_beta; nan(1,3)]; % trials x 3(2 regions + interactio)
                            All_NeurPairFF_pval = [All_NeurPairFF_pval; nan(1,3)];
                            All_NeurPairFF_Brain = [All_NeurPairFF_Brain; bregtmp];
                            All_NeurPairFF_NeurID = [All_NeurPairFF_NeurID; ...
                                [neurons_thisset(nn) neurons_thisset(nnn)]];
                            All_NeurPairFF_BirdID = [All_NeurPairFF_BirdID; ...
                                i];
                            All_NeurPairFF_ExptID = [All_NeurPairFF_ExptID; ...
                                ii];
                            All_NeurPairFF_Ntrials = [All_NeurPairFF_Ntrials; ...
                                length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract)];
                            All_NeurPairFF_motifnum = [All_NeurPairFF_motifnum; ...
                                mm];
                            All_NeurPairFF_setnum = [All_NeurPairFF_setnum; ...
                                iii];
                            continue
                        end
                        
                        X = NspksAll(:, [nn nnn]);
                        Y = FF;
                        
                        if tempfix ==1
                            if size(X,1)<5
                                continue
                            end
                        end
                        
                        
                        if handcodeinteraction==1
                            % ---- make new column by multiplying X (and then subtrctig mean)
                            % ---- since other way (using interactions) can
                            % have multiplication of negative numbers (if
                            % subtract mean first).
                            X = [X, X(:,1).* X(:,2)];
                            modelspec = 'linear';
                        else
                            modelspec = 'interactions';
                        end
                        
                        % --- subtract mean
                        X = X - repmat(mean(X,1), size(X,1), 1);
                        Y = FF - mean(FF);
                        
                        
                        mdl = fitlm(X, Y, modelspec);
                        
                        % ==== extract coeff, se, and p, from mdl
                        betas = mdl.Coefficients.Estimate(2:end);
                        pvals = mdl.Coefficients.pValue(2:end);
                        betasSE = mdl.Coefficients.SE(2:end);
                        
                        % ========= collect
                        All_NeurPairFF_beta = [All_NeurPairFF_beta; betas']; % trials x 3(2 regions + interactio)
                        All_NeurPairFF_pval = [All_NeurPairFF_pval; pvals'];
                        All_NeurPairFF_Brain = [All_NeurPairFF_Brain; bregtmp];
                        All_NeurPairFF_NeurID = [All_NeurPairFF_NeurID; ...
                            [neurons_thisset(nn) neurons_thisset(nnn)]];
                        All_NeurPairFF_BirdID = [All_NeurPairFF_BirdID; ...
                            i];
                        All_NeurPairFF_ExptID = [All_NeurPairFF_ExptID; ...
                            ii];
                        All_NeurPairFF_Ntrials = [All_NeurPairFF_Ntrials; ...
                            length(FF)];
                        All_NeurPairFF_motifnum = [All_NeurPairFF_motifnum; ...
                            mm];
                        All_NeurPairFF_setnum = [All_NeurPairFF_setnum; ...
                            iii];
                        
                        
                        % ================= PLOT EXAMPLES
                        if plotRaw==1
                            % -- region 1 vs FF
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            title(motifstr);
                            ylabel('FF');
                            xlabel(['n' num2str(neurons_thisset(nn)) '-' BrainregionsList{nn}]);
                            lt_regress(Y, X(:,1), 1);
                            lt_plot_zeroline;
                            lt_plot_zeroline_vert
                            
                            % -- region 2 vs. FF
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            title(motifstr);
                            ylabel('FF');
                            xlabel(['n' num2str(neurons_thisset(nnn)) '-' BrainregionsList{nnn}]);
                            lt_regress(Y, X(:,2), 1);
                            lt_plot_zeroline;
                            lt_plot_zeroline_vert
                            
                            
                            % -- region1*region2 vs. FF
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            title(motifstr);
                            ylabel('FF');
                            xlabel('interaction');
                            lt_regress(Y, X(:,1).*X(:,2), 1);
                            lt_plot_zeroline;
                            lt_plot_zeroline_vert
                            
                            %                             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            %                         title(motifstr);
                            %                         lt_plot_stem3(X(:,1), X(:,2), Y, 'r', 1);
                            %
                            %                         xlabel(['n' num2str(neurons_thisset(nn)) '-' BrainregionsList{nn}]);
                            %                         ylabel(['n' num2str(neurons_thisset(nnn)) '-' BrainregionsList{nnn}]);
                            %                         zlabel('FF');
                            %                         lt_plot_zeroline;
                            %                             lt_plot_zeroline_vert;
                        end
                    end
                end
                
                if plotRaw==1
                    pause
                    close all;
                end
                
                %% 2) EACH NEURON REGRESS WITH FF
                
                neurFF_brain = BrainregionsList;
                
                if isnan(FF(1))
                    All_NeurFF_Rho = [All_NeurFF_Rho nan(1,nneur)];
                    All_NeurFF_P = [All_NeurFF_P nan(1,nneur)];
                    All_NeurFF_Brain = [All_NeurFF_Brain neurFF_brain];
                    
                else
                    
                    [rho, p] = corr([NspksAll FF]);
                    neurFF_rho = rho(1:end-1,end); % corr between individual neurons and FF (in order)
                    neurFF_p = p(1:end-1, end);
                    
                    All_NeurFF_Rho = [All_NeurFF_Rho neurFF_rho'];
                    All_NeurFF_P = [All_NeurFF_P neurFF_p'];
                    All_NeurFF_Brain = [All_NeurFF_Brain neurFF_brain];
                    
                end
                
                %% 3) CORRELATION BETWEEN NEURONS
                for nn=1:nneur
                    for nnn=nn+1:nneur
                        
                        % ============================ SPIKE CORRELATIONS
                        rtmp = rho(nn, nnn);
                        ptmp = p(nn, nnn);
                        bregtmp = {BrainregionsList{nn}, BrainregionsList{nnn}};
                        
                        % ===== plot this as example?
                        if plotRaw_CC==1 & any(strcmp(bregtmp, 'LMAN')) & ...
                                any(strcmp(bregtmp, 'RA'))
                            plotThisCC = 1;
                        else
                            plotThisCC =0;
                        end
                        
                        
                        % ----- OUTPUT
                        All_NeurNeur_Rho = [All_NeurNeur_Rho; rtmp'];
                        All_NeurNeur_P = [All_NeurNeur_P; ptmp'];
                        All_NeurNeur_Brain = [All_NeurNeur_Brain; bregtmp];
                        
                        
                        % ============================ CROSS CORRELATIONS
                        dattmp1 = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(nn);
                        dattmp2 = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(nnn);
                        
                        clustnum1 = SummaryStruct.birds(i).neurons(neurons_thisset(nn)).clustnum;
                        clustnum2 = SummaryStruct.birds(i).neurons(neurons_thisset(nnn)).clustnum;
                        
                        % -------- extract FR
                        dattmp1.SegmentsExtract = lt_neural_SmoothFR(dattmp1.SegmentsExtract, ...
                            clustnum1, '', '', '', segextract_for_trialdur);
                        dattmp2.SegmentsExtract = lt_neural_SmoothFR(dattmp2.SegmentsExtract, ...
                            clustnum2, '', '', '', segextract_for_trialdur);
                        
                        
                        plotcol1 = plotCols{strcmp(plotColsRegions, bregtmp{1})};
                        plotcol2 = plotCols{strcmp(plotColsRegions, bregtmp{2})};
                        
                        % ============================ CROSS CORRELATION
                        ntrials = length(dattmp1.SegmentsExtract);
                        CCall = [];
                        CCall_shift = [];
                        FRall1 = []; % first region
                        FRall2 = []; % second region
                        
                        if plotThisCC==1
                            lt_figure; hold on;
                            ypos = 1;
                            
                            lt_subplot(8, 2, 1:14); hold on;
                            title([bregtmp{1} '-' bregtmp{2} '[' motifstr ']-set' ...
                                num2str(iii) '-n' num2str(neurons_thisset(nn)) ',' ...
                                num2str(neurons_thisset(nnn))]);
                        end
                        
                        for t = 1:ntrials
                            
                            fr1 = dattmp1.SegmentsExtract(t).FRsmooth_rate_CommonTrialDur;
                            fr2 = dattmp2.SegmentsExtract(t).FRsmooth_rate_CommonTrialDur;
                            frx = dattmp1.SegmentsExtract(t).FRsmooth_xbin_CommonTrialDur;
                            
                            % --- limit to data in "premotor" window
                            windtmp = motifpredur + xcov_dattotake;
                            fr1_short = fr1(frx>=windtmp(1) & frx<windtmp(2));
                            fr2_short = fr2(frx>=windtmp(1) & frx<windtmp(2));
                            
                            
                            % ------------------- calculate
                            [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
                            %                             [cc, lags] = xcov(fr1, fr2, xcovwindmax/0.001, 'coeff');
                            
                            CCall = [CCall; cc'];
                            FRall1 = [FRall1; fr1_short'];
                            FRall2 = [FRall2; fr2_short'];
                            
                            % ############################ PLOT RAW
                            if plotThisCC ==1
                                % ----- plot random subset of trials, plot
                                % each channel overlaied (rasters + FR)
                                if rand < 1/(ntrials/15)
                                    lt_subplot(8, 2, 1:14); hold on;
                                    
                                    % 1) =================== raster
                                    spktimes1 = dattmp1.SegmentsExtract(t).spk_Times(...
                                        dattmp1.SegmentsExtract(t).spk_Clust == clustnum1);
                                    spktimes2 = dattmp2.SegmentsExtract(t).spk_Times(...
                                        dattmp2.SegmentsExtract(t).spk_Clust == clustnum2);
                                    
                                    
                                    lt_neural_PLOT_rasterline(spktimes1, ypos, plotcol1, 0);
                                    lt_neural_PLOT_rasterline(spktimes2, ypos+1, plotcol2, 0);
                                    
                                    
                                    % 2) ===================== smoothed FR
                                    %                                  lt_subplot(1,2,2); % smoothed FR
                                    x = (1:length(fr1))*0.001;
                                    plot(x, (ypos-0.4)+fr1./(max(fr1)), '-k', 'LineWidth', 2);
                                    plot(x, (ypos+0.4)+fr2./(max(fr2)), '-k', 'LineWidth', 2);
                                    
                                    % 3) =============== syl contour
                                    
                                    lt_plot(segextract_for_trialdur(t).sylOnTimes_RelDataOnset, ypos-0.5, {'Color', 'g'});
                                    lt_plot(segextract_for_trialdur(t).sylOffTimes_RelDataOnset, ypos-0.5, {'Color', 'r'});
                                    
                                    
                                    % ------------
                                    ypos = ypos+3;
                                end
                                
                                
                            end
                            
                            
                            % #################################
                            
                            
                            % ----------------- shift control
                            if t<ntrials
                                fr2 = dattmp2.SegmentsExtract(t+1).FRsmooth_rate_CommonTrialDur;
                            else
                                fr2 = dattmp2.SegmentsExtract(1).FRsmooth_rate_CommonTrialDur;
                            end
                            fr2_short = fr2(frx>=windtmp(1) & frx<windtmp(2));
                            
                            
                            [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
                            CCall_shift = [CCall_shift; cc'];
                            
                        end
                        
                        if plotThisCC ==1
                            
                            % ---- line showing premotor window
                            lt_subplot(8,2,1:14); hold on;
                            axis tight;
                            line([windtmp(1) windtmp(1)], ylim, 'Color', 'k');
                            line([windtmp(2) windtmp(2)], ylim, 'Color', 'k');
                            
                            % ----
                            lt_subplot(8, 2, 15); hold on;
                            title('cross corr');
                            
                            shadedErrorBar(lags*0.001, mean(CCall_shift,1), lt_sem(CCall_shift), {'Color', [0.7 0.7 0.7]}, 1);
                            shadedErrorBar(lags*0.001, mean(CCall,1), lt_sem(CCall), {'Color', 'k'}, 1);
                            
                            %                             lt_plot(lags*0.001, mean(CCall_shift,1), {'Errors', lt_sem(CCall_shift), 'Color', [0.7 0.7 0.7]});
                            %                             lt_plot(lags*0.001, mean(CCall,1), {'Errors', lt_sem(CCall), 'Color', 'k'});
                            
                            
                            pause
                            close all;
                            
                        end
                        
                        
                        % ----- OUTPUT
                        
                        All_NeurNeur_xcorr_shift = [All_NeurNeur_xcorr_shift; CCall_shift];
                        All_NeurNeur_xcorr = [All_NeurNeur_xcorr; CCall];
                        All_NeurNeur_xcorr_lags = [All_NeurNeur_xcorr_lags; lags];
                        All_NeurNeur_FRreg1 = [All_NeurNeur_FRreg1; FRall1];
                        All_NeurNeur_FRreg2 = [All_NeurNeur_FRreg2; FRall2];
                        
                        
                    end
                end
                
                
                
                if (0)
                    
                    lt_figure; hold on;
                    plot(NspksAll(:,1).*NspksAll(:,2), FF, 'ok');
                    
                end
            end
        end
    end
end

%%
inds = All_NeurPairFF_BirdID==1 & All_NeurPairFF_motifnum==1;
All_NeurPairFF_NeurID(inds, :)

%%
tmp = strcmp(All_NeurPairFF_Brain, All_NeurNeur_Brain);

assert(all(tmp(:)==1)==1, 'asdfasd');

%% ############################## make brain regions alphabetical order
% reason: so region 1 is always same actual region, for a given pair

for j=1:size(All_NeurPairFF_Brain, 1)
    
    [~, indtmp] = sort(All_NeurPairFF_Brain(j,:));
    
    % -- output
    All_NeurPairFF_Brain(j,:) = All_NeurPairFF_Brain(j, indtmp);
    All_NeurPairFF_beta(j, [1 2]) = All_NeurPairFF_beta(j, indtmp);
    All_NeurPairFF_pval(j, [1 2]) = All_NeurPairFF_pval(j, indtmp);
    
    % --
    clear All_NeurNeur_Brain; % replaced by All_NeurPairFF_Brain
    if indtmp(2)>indtmp(1)
        % -- flip xcorr lr
        All_NeurNeur_xcorr{j} = fliplr(All_NeurNeur_xcorr{j});
        All_NeurNeur_xcorr_shift{j} = fliplr(All_NeurNeur_xcorr_shift{j});
        
        tmp = All_NeurNeur_FRreg1{j};
        All_NeurNeur_FRreg1{j} = All_NeurNeur_FRreg2{j};
        All_NeurNeur_FRreg2{j} = tmp;
        
    end
end

%%
BrainregionsList = unique(All_NeurFF_Brain);


%% ############################## PLOTS

%% =============== PLOT ALL MEAN CROSS CORRELATIONS

assert(all(diff(All_NeurPairFF_NeurID,1 ,2)>0), 'asfasdf');

maxExpts = max(All_NeurPairFF_ExptID);
maxNeur = max(All_NeurPairFF_NeurID(:));
maxBirds = max(All_NeurPairFF_BirdID);
maxSet = max(All_NeurPairFF_setnum);

for i=1:maxBirds
    
    FR1all = [];
    FR2all = [];    
    CCall = [];
    CCshiftall = [];
    NeurPairIDall = [];
    BregionPairAll = {};
    MotifNumAll = [];
    npairID = 0;
    
    for ii=1:maxExpts
        
        figstart = 0;
        
            for nn = 1:maxNeur
                for nnn = nn+1:maxNeur
                    
                    inds = All_NeurPairFF_BirdID==i & All_NeurPairFF_ExptID==ii ...
                        & All_NeurPairFF_NeurID(:,1)==nn & All_NeurPairFF_NeurID(:,2)==nnn;
                    
                    if ~any(inds)
                        continue
                    end
                    
                    assert(length(unique(All_NeurPairFF_setnum(inds)))==1, 'neuron pair in multiple sets. would be redundant to count multiple times. need to write code to combine i.e. add trials');
                    
                    % -- increment neur pair ID
                    npairID = npairID+1;
                   
                    if figstart ==0
                        % --- initiate figure
                        figcount=1;
                        subplotrows=5;
                        subplotcols=sum(inds); % number of motifs
                        fignums_alreadyused=[];
                        hfigs=[];
                        hsplots = [];
                        
                        figstart=1;
                    end
                    
                    % -- sanity check
                    if any(inds) & (0)
                        disp('------------');
                        disp([num2str(nn) '-' num2str(nnn)]);
                        disp([All_NeurPairFF_motifnum(inds) All_NeurPairFF_setnum(inds)]);
                        
                    end
                    
                    
                    % ======== extract all data for this neuron pair
                    bregion1 = unique(All_NeurPairFF_Brain(inds,1));
                    bregion2 = unique(All_NeurPairFF_Brain(inds,2));
                    bregion1 = bregion1{1};
                    bregion2 = bregion2{1};
                    bregionpair = [bregion1 '-' bregion2];
                    
                    xlags = mean(All_NeurNeur_xcorr_lags(inds, :),1);
                    motifsThis = All_NeurPairFF_motifnum(inds);
                    ccThis = All_NeurNeur_xcorr(inds);
                    ccShiftThis = All_NeurNeur_xcorr_shift(inds);
                    
                    frthis1 = All_NeurNeur_FRreg1(inds);
                    frthis2 = All_NeurNeur_FRreg2(inds);
                    
                    % ======= plot
                    for mm = 1:length(motifsThis)
                        motifnum = motifsThis(mm);
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title([MOTIFSTATS_pop.birds(i).params.motif_regexpr_str{motifnum}]);
                        xlabel([bregion1 '-' bregion2]);
                        if mm==1
                            ylabel(['nn' num2str(nn) '-' num2str(nnn)]);
                        end
                        
                        % ------- plot
                        % - shift
                        ymat = ccShiftThis{mm};
                        %                     ymat = cell2mat(cellfun(@nanmean, ccShiftThis, 'UniformOutput', 0));
                        ymean = nanmean(ymat);
                        ysem = lt_sem(ymat);
                        shadedErrorBar(xlags*0.001, ymean, ysem, {'Color', [0.7 0.4 0.4]}, 1);
                        
                        %                     ymat = cell2mat(cellfun(@nanmean, ccThis, 'UniformOutput', 0));
                        ymat = ccThis{mm};
                        ymean = nanmean(ymat);
                        ysem = lt_sem(ymat);
                        shadedErrorBar(xlags*0.001, ymean, ysem, {'Color', 'k'}, 1);
                        
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                        
                                              
                        
                        % =============== collect for summary plot
                        CCall = [CCall; nanmean(ccThis{mm})];
                        CCshiftall = [CCshiftall; nanmean(ccShiftThis{mm})];
                        
                        FR1all = [FR1all; mean(frthis1{mm})];
                        FR2all = [FR2all; mean(frthis2{mm})];
                        
                        NeurPairIDall = [NeurPairIDall; npairID];
                        BregionPairAll = [BregionPairAll; bregionpair];
                        MotifNumAll = [MotifNumAll; motifnum];

                        
                        
                    end
                    
                end
        end
        linkaxes(hsplots, 'xy');
        lt_subtitle(['b' num2str(i) ',ex' num2str(ii)]);
        
    end
    
    % ################################## Summary plot for this bird [CC and
    % shift
    figcount=1;
    subplotrows=3;
    subplotcols=max(MotifNumAll); % number of motifs
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    bregionpairlist = {'LMAN-LMAN', 'LMAN-RA', 'RA-RA'};
    
    for bb = bregionpairlist
        
        % --- one plot for each motif
        for mm = 1:max(MotifNumAll)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            
            inds = MotifNumAll==mm & strcmp(BregionPairAll, bb);
            
            if mm==1
                ylabel(bb);
            end
            title(MOTIFSTATS_pop.birds(i).params.motif_regexpr_str{mm})
            % for this motif and brain region pair
            ccthis = CCall(inds, :);
            ccshiftthis = CCshiftall(inds, :);
            
            neurpairids = NeurPairIDall(inds);
            
            shadedErrorBar(xlags*0.001, nanmean(ccthis,1), lt_sem(ccthis), {'Color', 'k'}, 1);
            shadedErrorBar(xlags*0.001, nanmean(ccshiftthis,1), lt_sem(ccshiftthis), {'Color', [0.8 0.5 0.5]}, 1);
            
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
        end
    end
    
    linkaxes(hsplots, 'xy');
    xlim([-0.075 0.075]);
    ylim([-0.3 0.3])
    
    
    
        % ################################## Summary plot for this bird [CC
        % minus shift]
    figcount=1;
    subplotrows=3;
    subplotcols=max(MotifNumAll); % number of motifs
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    bregionpairlist = {'LMAN-LMAN', 'LMAN-RA', 'RA-RA'};
    
    for bb = bregionpairlist
        
        % --- one plot for each motif
        for mm = 1:max(MotifNumAll)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            
            inds = MotifNumAll==mm & strcmp(BregionPairAll, bb);
            
            if mm==1
                ylabel(bb);
            end
            title(MOTIFSTATS_pop.birds(i).params.motif_regexpr_str{mm})
            % for this motif and brain region pair
            ccthis = CCall(inds, :);
            ccshiftthis = CCshiftall(inds, :);
            
            ccminus = ccthis - ccshiftthis;
            
            plot(xlags*0.001, ccminus', 'Color', [0.7 0.3 0.3]);
%             lt_plot(xlags*0.001, nanmean(ccminus,1), {'Errors', lt_sem(ccminus), 'Color', 'k'});
            shadedErrorBar(xlags*0.001, nanmean(ccminus,1), lt_sem(ccminus), {'Color', 'b'}, 1);

            
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
        end
    end
    
    linkaxes(hsplots, 'xy');
    xlim([-0.075 0.075]);
    ylim([-0.25 0.25])
    
    % ############################################## FIRING RATE
    figcount=1;
    subplotrows=3;
    subplotcols=max(MotifNumAll); % number of motifs
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    bregionpairlist = {'LMAN-LMAN', 'LMAN-RA', 'RA-RA'};
    
    for bb = bregionpairlist
        
        % --- one plot for each motif
        for mm = 1:max(MotifNumAll)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            
            inds = MotifNumAll==mm & strcmp(BregionPairAll, bb);
            
            if mm==1
                ylabel([bb 'bk-rd']);
            end
            title(MOTIFSTATS_pop.birds(i).params.motif_regexpr_str{mm})
            % for this motif and brain region pair
            fr1 = FR1all(inds,:);
            fr2 = FR2all(inds,:);
            
            x = 0.001*(1:size(fr1,2));
            shadedErrorBar(x, mean(fr1,1), lt_sem(fr1), {'Color', 'k'}, 1);
            shadedErrorBar(x, mean(fr2,1), lt_sem(fr2), {'Color', 'r'}, 1);
            
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
        end
    end
    
    linkaxes(hsplots, 'xy');
    axis tight;
%     xlim([-0.075 0.075]);
%     ylim([-0.25 0.25])
% 
%     
end


%% ============== SINGLE NEURON VS. FF [DISTRIBUTION PER REGION]

for i=1:length(BrainregionsList)
    lt_figure; hold on;
    ct = 1;
    bregion = BrainregionsList{i};
    
    inds = strcmp(All_NeurFF_Brain, bregion) & ~isnan(All_NeurFF_Rho);
    
    % -------------- 1) histogram of rho
    lt_subplot(3,2,ct); hold on;
    ct = ct+1;
    title(bregion);
    xlabel('rho (all neurons)');
    
    rho = All_NeurFF_Rho(inds);
    p = All_NeurFF_P(inds);
    
    lt_plot_histogram(rho);
    xlim([-1 1]);
    
    % --------------- 2) rho vs. p value
    lt_subplot(3,2,ct); ct= ct+1; hold on;
    ylabel('log10(p)');
    xlabel('rho');
    
    plot(rho, log10(p), 'or');
    line(xlim, [log10(0.05) log10(0.05)]);
    
    % -------------- 3) count of significant negative and positive)
    lt_subplot(3,2,ct); ct= ct+1; hold on;
    xlabel('p threshold');
    ylabel('% pass');
    title(bregion);
    
    propsig_pos = [];
    propsig_neg = [];
    threshlist = [0.05 0.01 0.001];
    for thresh =threshlist;
        % -- positive
        propsig_pos = [propsig_pos sum(p<thresh & rho>0)/length(p)];
        % -- neg
        propsig_neg = [propsig_neg sum(p<thresh & rho<0)/length(p)];
        
    end
    % -- pos
    lt_plot_bar([1:length(threshlist)]-0.2, propsig_pos, {'BarWidth', 0.2, 'Color', 'b'});
    % -- neg
    lt_plot_bar([1:length(threshlist)]+0.2, propsig_neg, {'BarWidth', 0.2, 'Color', 'r'});
    
    xlabel(num2str(threshlist));
    lt_plot_annotation(1, ['N=' num2str(length(p))], 'r');
    ylim([0 1]);
end

%% ============= NEURON PAIRS VS. FF
numrows = 3; numcols = 2;

for i = 1:length(BrainregionsList)
    bregion1 = BrainregionsList{i};
    
    % ##################################### THIS REGION PAIRED WITH ITSELF
    lt_figure; hold on;
    ct = 1;
    
    inds = strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
        & strcmp(All_NeurPairFF_Brain(:,2), bregion1) ...
        & ~isnan(All_NeurPairFF_beta(:,1));
    
    % ================================= 1) histogram of beta
    lt_subplot(numrows,numcols,ct); hold on;
    ct = ct+1;
    title([bregion1 '-' bregion1]);
    xlabel('betas');
    
    % --- first region
    beta = All_NeurPairFF_beta(inds, 1);
    plotcol = 'r';
    lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
    lt_plot_text(max(beta), 2, 'region1', plotcol);
    
    % --- second region
    beta = All_NeurPairFF_beta(inds, 2);
    plotcol = 'b';
    lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
    lt_plot_text(max(beta), 4, 'region2', plotcol);
    
    % --- interaction
    beta = All_NeurPairFF_beta(inds, 3);
    plotcol = 'k';
    lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
    lt_plot_text(max(beta), 6, 'interaction', plotcol);
    
    
    % ================================ 3D scatter plot
    lt_subplot(numrows,numcols,ct); hold on;
    ct = ct+1;
    xlabel('region1', 'Color', 'r');
    ylabel('region2', 'Color', 'b');
    zlabel('interaction', 'Color', 'k');
    
    % -- all
    beta1 = All_NeurPairFF_beta(inds, 1);
    beta2 = All_NeurPairFF_beta(inds, 2);
    beta3 = All_NeurPairFF_beta(inds, 3);
    %         plot3(beta1, beta2, beta3, 'ok');
    lt_plot_stem3(beta1, beta2, beta3, 'k', 1);
    
    
    
    % ================================ 2) beta vs. p value
    lt_subplot(numrows,numcols,ct); hold on;
    ct = ct+1;
    ylabel('log10(p)');
    xlabel('beta');
    
    
    % --- first region
    beta = All_NeurPairFF_beta(inds, 1);
    p = All_NeurPairFF_pval(inds, 1);
    plotcol = 'r';
    
    plot(beta, log10(p), 'o', 'Color', plotcol);
    
    % --- second region
    beta = All_NeurPairFF_beta(inds, 2);
    p = All_NeurPairFF_pval(inds, 2);
    plotcol = 'b';
    
    plot(beta, log10(p), '^', 'Color', plotcol);
    
    % --- third region
    beta = All_NeurPairFF_beta(inds, 3);
    p = All_NeurPairFF_pval(inds, 3);
    plotcol = 'k';
    
    plot(beta, log10(p), 's', 'Color', plotcol);
    
    
    line(xlim, [log10(0.05) log10(0.05)]);
    lt_plot_zeroline_vert;
    
    % ======================== 3) count of significant negative and positive)
    % ------ iterate over region 1, region 2, interaction
    for j=1:3;
        lt_subplot(numrows,numcols,ct); hold on; ct = ct+1;
        xlabel('p threshold');
        ylabel('% pass');
        if j==3
            title('interaction');
        else
            title(['region' num2str(j)])
        end
        
        p = All_NeurPairFF_pval(inds, j);
        beta = All_NeurPairFF_beta(inds, j);
        
        propsig_pos = [];
        propsig_neg = [];
        threshlist = [0.05 0.01 0.001];
        for thresh =threshlist;
            % -- positive
            propsig_pos = [propsig_pos sum(p<thresh & beta>0)/length(p)];
            % -- neg
            propsig_neg = [propsig_neg sum(p<thresh & beta<0)/length(p)];
        end
        
        % -- pos
        lt_plot_bar([1:length(threshlist)]-0.2, propsig_pos, {'BarWidth', 0.2, 'Color', 'b'});
        % -- neg
        lt_plot_bar([1:length(threshlist)]+0.2, propsig_neg, {'BarWidth', 0.2, 'Color', 'r'});
        
        xlabel(num2str(threshlist));
        lt_plot_annotation(1, ['N=' num2str(length(p))], 'r');
        ylim([0 1]);
    end
    
    
    % ############################################### BETWEEN BRAIN REGIONS
    for ii=i+1:length(BrainregionsList)
        bregion2 = BrainregionsList{ii};
        lt_figure; hold on;
        ct = 1;
        
        inds = ((strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
            & strcmp(All_NeurPairFF_Brain(:,2), bregion2)) ...
            | ...
            (strcmp(All_NeurPairFF_Brain(:,1), bregion2) ...
            & strcmp(All_NeurPairFF_Brain(:,2), bregion1))) ...
            &  ~isnan(All_NeurPairFF_beta(:,1)); % since bould be 1-2 or 2-1 (diff order)
        
        
        % ================================= 1) histogram of beta
        lt_subplot(numrows,numcols,ct); hold on;
        ct = ct+1;
        title([bregion1 '-' bregion2]);
        xlabel('betas');
        
        % --- first region
        beta = All_NeurPairFF_beta(inds, 1);
        plotcol = 'r';
        lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
        lt_plot_text(max(beta), 2, 'region1', plotcol);
        
        % --- second region
        beta = All_NeurPairFF_beta(inds, 2);
        plotcol = 'b';
        lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
        lt_plot_text(max(beta), 4, 'region2', plotcol);
        
        % --- interaction
        beta = All_NeurPairFF_beta(inds, 3);
        plotcol = 'k';
        lt_plot_histogram(beta, '', 1, 0, '', 1, plotcol);
        lt_plot_text(max(beta), 6, 'interaction', plotcol);
        
        
        % ================================ 3D scatter plot
        lt_subplot(numrows,numcols,ct); hold on;
        ct = ct+1;
        xlabel('region1', 'Color', 'r');
        ylabel('region2', 'Color', 'b');
        zlabel('interaction', 'Color', 'k');
        
        % -- all
        beta1 = All_NeurPairFF_beta(inds, 1);
        beta2 = All_NeurPairFF_beta(inds, 2);
        beta3 = All_NeurPairFF_beta(inds, 3);
        %         plot3(beta1, beta2, beta3, 'ok');
        lt_plot_stem3(beta1, beta2, beta3, 'k', 1);
        
        % ================================ 2) beta vs. p value
        lt_subplot(numrows,numcols,ct); hold on;
        ct = ct+1;
        ylabel('log10(p)');
        xlabel('beta');
        
        
        % --- first region
        beta = All_NeurPairFF_beta(inds, 1);
        p = All_NeurPairFF_pval(inds, 1);
        plotcol = 'r';
        
        plot(beta, log10(p), 'o', 'Color', plotcol);
        
        % --- second region
        beta = All_NeurPairFF_beta(inds, 2);
        p = All_NeurPairFF_pval(inds, 2);
        plotcol = 'b';
        
        plot(beta, log10(p), '^', 'Color', plotcol);
        
        % --- third region
        beta = All_NeurPairFF_beta(inds, 3);
        p = All_NeurPairFF_pval(inds, 3);
        plotcol = 'k';
        
        plot(beta, log10(p), 's', 'Color', plotcol);
        
        
        line(xlim, [log10(0.05) log10(0.05)]);
        lt_plot_zeroline_vert;
        
        % ======================== 3) count of significant negative and positive)
        % ------ iterate over region 1, region 2, interaction
        for j=1:3;
            lt_subplot(numrows,numcols,ct); hold on; ct = ct+1;
            xlabel('p threshold');
            ylabel('% pass');
            if j==3
                title('interaction');
            else
                title(['region' num2str(j)])
            end
            
            p = All_NeurPairFF_pval(inds, j);
            beta = All_NeurPairFF_beta(inds, j);
            
            propsig_pos = [];
            propsig_neg = [];
            threshlist = [0.05 0.01 0.001];
            for thresh =threshlist;
                % -- positive
                propsig_pos = [propsig_pos sum(p<thresh & beta>0)/length(p)];
                % -- neg
                propsig_neg = [propsig_neg sum(p<thresh & beta<0)/length(p)];
            end
            
            % -- pos
            lt_plot_bar([1:length(threshlist)]-0.2, propsig_pos, {'BarWidth', 0.2, 'Color', 'b'});
            % -- neg
            lt_plot_bar([1:length(threshlist)]+0.2, propsig_neg, {'BarWidth', 0.2, 'Color', 'r'});
            
            xlabel(num2str(threshlist));
            lt_plot_annotation(1, ['N=' num2str(length(p))], 'r');
            ylim([0 1]);
        end
        
        
    end
end

%% ============= NEURON-NEURON PAIRS

for i = 1:length(BrainregionsList)
    bregion1 = BrainregionsList{i};
    % ============= THIS REGION PAIRED WITH ITSELF
    lt_figure; hold on;
    ct = 1;
    
    inds = strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
        & strcmp(All_NeurPairFF_Brain(:,2), bregion1);
    
    % ------------------------------------------------ PLOTS
    % -------------- 1) histogram of rho
    lt_subplot(2,2,ct); hold on;
    ct = ct+1;
    title([bregion1 '-' bregion1]);
    xlabel('rho (all neurons)');
    
    rho = All_NeurNeur_Rho(inds);
    p = All_NeurNeur_P(inds);
    
    lt_plot_histogram(rho);
    xlim([-1 1]);
    
    % --------------- 2) rho vs. p value
    lt_subplot(2,2,ct); ct= ct+1; hold on;
    ylabel('log10(p)');
    xlabel('rho');
    
    plot(rho, log10(p), 'or');
    line(xlim, [log10(0.05) log10(0.05)]);
    
    
    % -------------- 3) count of significant negative and positive)
    lt_subplot(2,2,ct); ct= ct+1; hold on;
    xlabel('p threshold');
    ylabel('% pass');
    
    propsig_pos = [];
    propsig_neg = [];
    threshlist = [0.05 0.01 0.001];
    for thresh =threshlist;
        % -- positive
        propsig_pos = [propsig_pos sum(p<thresh & rho>0)/length(p)];
        % -- neg
        propsig_neg = [propsig_neg sum(p<thresh & rho<0)/length(p)];
    end
    
    % -- pos
    lt_plot_bar([1:length(threshlist)]-0.2, propsig_pos, {'BarWidth', 0.2, 'Color', 'b'});
    % -- neg
    lt_plot_bar([1:length(threshlist)]+0.2, propsig_neg, {'BarWidth', 0.2, 'Color', 'r'});
    
    xlabel(num2str(threshlist));
    lt_plot_annotation(1, ['N=' num2str(length(p))], 'r');
    ylim([0 1]);
    
    % ========================= THIS REGION PAIRED WITH OTHERS
    for ii=i+1:length(BrainregionsList)
        lt_figure; hold on;
        ct = 1;
        
        bregion2 = BrainregionsList{ii};
        
        inds = (strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
            & strcmp(All_NeurPairFF_Brain(:,2), bregion2)) ...
            | ...
            (strcmp(All_NeurPairFF_Brain(:,1), bregion2) ...
            & strcmp(All_NeurPairFF_Brain(:,2), bregion1)); % since bould be 1-2 or 2-1 (diff order)
        
        % ------------------------------------------------ PLOTS
        % -------------- 1) histogram of rho
        lt_subplot(2,2,ct); hold on;
        ct = ct+1;
        title([bregion1 '-' bregion2]);
        xlabel('rho (all neurons)');
        
        rho = All_NeurNeur_Rho(inds);
        p = All_NeurNeur_P(inds);
        
        lt_plot_histogram(rho);
        xlim([-1 1]);
        
        % --------------- 2) rho vs. p value
        lt_subplot(2,2,ct); ct= ct+1; hold on;
        ylabel('log10(p)');
        xlabel('rho');
        
        plot(rho, log10(p), 'or');
        line(xlim, [log10(0.05) log10(0.05)]);
        
        
        % -------------- 3) count of significant negative and positive)
        lt_subplot(2,2,ct); ct= ct+1; hold on;
        xlabel('p threshold');
        ylabel('% pass');
        
        propsig_pos = [];
        propsig_neg = [];
        threshlist = [0.05 0.01 0.001];
        for thresh =threshlist;
            % -- positive
            propsig_pos = [propsig_pos sum(p<thresh & rho>0)/length(p)];
            % -- neg
            propsig_neg = [propsig_neg sum(p<thresh & rho<0)/length(p)];
        end
        
        % -- pos
        lt_plot_bar([1:length(threshlist)]-0.2, propsig_pos, {'BarWidth', 0.2, 'Color', 'b'});
        % -- neg
        lt_plot_bar([1:length(threshlist)]+0.2, propsig_neg, {'BarWidth', 0.2, 'Color', 'r'});
        
        xlabel(num2str(threshlist));
        lt_plot_annotation(1, ['N=' num2str(length(p))], 'r');
        ylim([0 1]);
        
        
    end
end


%% ================== CROSS CORRELATION BETWEEN BRAIN REGIONS


for i = 1:length(BrainregionsList)
    bregion1 = BrainregionsList{i};
    % ============= THIS REGION PAIRED WITH ITSELF
    lt_figure; hold on;
    ct = 1;
    
    inds = strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
        & strcmp(All_NeurPairFF_Brain(:,2), bregion1);
    
    % ------------------------------------------------ PLOTS
    % -------------- 1) histogram of rho
    lt_subplot(2,2,ct); hold on;
    ct = ct+1;
    title('xcorr of FR');
    ylabel('mean (across trials) of xcorr');
    xlabel([bregion1 '<--->' bregion1]);
    %     xlabel('cc (all neurons)');
    
    ymean = cellfun(@nanmean, All_NeurNeur_xcorr(inds), 'UniformOutput', 0); % mean across trials
    ymean = cell2mat(ymean);
    
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
    
    
    % -------------- 1) histogram of rho [SHIFTED]
    lt_subplot(2,2,ct); hold on;
    ct = ct+1;
    title('xcorr of FR (SHIFT)');
    ylabel('mean (across trials) of xcorr');
    xlabel([bregion1 '<--->' bregion1]);
    %     xlabel('cc (all neurons)');
    
    ymean = cellfun(@nanmean, All_NeurNeur_xcorr_shift(inds), 'UniformOutput', 0); % mean across trials
    ymean = cell2mat(ymean);
    
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
    
    
    % ----------------- SUBTRACT SHIFT FROM ACTUAL
    lt_subplot(2,2,ct); hold on;
    ct = ct+1;
    title('xcorr of FR (dat minus shift)');
    ylabel('mean (across trials) of xcorr');
    xlabel([bregion1 '<--->' bregion1]);
    %     xlabel('cc (all neurons)');
    
    ymeanDAT = cellfun(@nanmean, All_NeurNeur_xcorr(inds), 'UniformOutput', 0); % mean across trials
    ymeanDAT = cell2mat(ymeanDAT);
    
    ymeanSHIFT = cellfun(@nanmean, All_NeurNeur_xcorr_shift(inds), 'UniformOutput', 0); % mean across trials
    ymeanSHIFT = cell2mat(ymeanSHIFT);
    
    ymean = ymeanDAT - ymeanSHIFT;
    
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
    plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
    
    
    % ========================= THIS REGION PAIRED WITH OTHERS
    for ii=i+1:length(BrainregionsList)
        lt_figure; hold on;
        ct = 1;
        
        bregion2 = BrainregionsList{ii};
        
        % --- put brain region names in order
        br = sort({bregion1, bregion2});
        bregion1 = br{1};
        bregion2 = br{2};
        
        % ---
        inds = (strcmp(All_NeurPairFF_Brain(:,1), bregion1) ...
            & strcmp(All_NeurPairFF_Brain(:,2), bregion2));
        
        % ------------------------------------------------ PLOTS
        % -------------- 1) histogram of rho
        lt_subplot(2,2,ct); hold on;
        ct = ct+1;
        title('xcorr of FR');
        ylabel('mean (across trials) of xcorr');
        xlabel([bregion1 '<--->' bregion2]);
        %     xlabel('cc (all neurons)');
        
        ymean = cellfun(@nanmean, All_NeurNeur_xcorr(inds), 'UniformOutput', 0); % mean across trials
        ymean = cell2mat(ymean);
        
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
        
        
        % -------------- 1) histogram of rho [SHIFTED]
        lt_subplot(2,2,ct); hold on;
        ct = ct+1;
        title('xcorr of FR (SHIFT)');
        ylabel('mean (across trials) of xcorr');
        xlabel([bregion1 '<--->' bregion2]);
        %     xlabel('cc (all neurons)');
        
        ymean = cellfun(@nanmean, All_NeurNeur_xcorr_shift(inds), 'UniformOutput', 0); % mean across trials
        ymean = cell2mat(ymean);
        
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
        
        
        % ----------------- SUBTRACT SHIFT FROM ACTUAL
        lt_subplot(2,2,ct); hold on;
        ct = ct+1;
        title('xcorr of FR (dat minus shift)');
        ylabel('mean (across trials) of xcorr');
        xlabel([bregion1 '<--->' bregion2]);
        %     xlabel('cc (all neurons)');
        
        ymeanDAT = cellfun(@nanmean, All_NeurNeur_xcorr(inds), 'UniformOutput', 0); % mean across trials
        ymeanDAT = cell2mat(ymeanDAT);
        
        ymeanSHIFT = cellfun(@nanmean, All_NeurNeur_xcorr_shift(inds), 'UniformOutput', 0); % mean across trials
        ymeanSHIFT = cell2mat(ymeanSHIFT);
        
        ymean = ymeanDAT - ymeanSHIFT;
        
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, ymean', 'k')
        plot(All_NeurNeur_xcorr_lags(1,:)*0.001, mean(ymean), 'r')
        
    end
end
