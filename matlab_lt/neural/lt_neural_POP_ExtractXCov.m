function MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct)
%% extracts pairwise cross correlations. sticks back into structure

%%
onlySetsWithRALMANpaired = 0;
plotRaw_CC = 0;

plotCols = lt_make_plot_colors(3, 0, 0);
plotColsRegions = {'LMAN', 'RA', 'X'};


% ========== xcorr analysis of frate
xcovwindmax = 0.05;
xcov_dattotake = [-0.1 0.05]; % rel syl onset

%%

NumBirds = length(MOTIFSTATS_pop.birds);

for i=1:NumBirds
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);

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
                
                
%                 %% 1) PERFORM MULTILINEAR REGRESSION [ALL PAIRS]
%                 for nn=1:nneur
%                     for nnn=nn+1:nneur
%                         
%                         % --- brain regions
%                         bregtmp = {BrainregionsList{nn}, BrainregionsList{nnn}};
%                         
%                         
%                         
%                         if isnan(FF(1))
%                             % then skip this neuron pair
%                             All_NeurPairFF_beta = [All_NeurPairFF_beta; nan(1,3)]; % trials x 3(2 regions + interactio)
%                             All_NeurPairFF_pval = [All_NeurPairFF_pval; nan(1,3)];
%                             All_NeurPairFF_Brain = [All_NeurPairFF_Brain; bregtmp];
%                             All_NeurPairFF_NeurID = [All_NeurPairFF_NeurID; ...
%                                 [neurons_thisset(nn) neurons_thisset(nnn)]];
%                             All_NeurPairFF_BirdID = [All_NeurPairFF_BirdID; ...
%                                 i];
%                             All_NeurPairFF_ExptID = [All_NeurPairFF_ExptID; ...
%                                 ii];
%                             All_NeurPairFF_Ntrials = [All_NeurPairFF_Ntrials; ...
%                                 length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract)];
%                             All_NeurPairFF_motifnum = [All_NeurPairFF_motifnum; ...
%                                 mm];
%                             All_NeurPairFF_setnum = [All_NeurPairFF_setnum; ...
%                                 iii];
%                             continue
%                         end
%                         
%                         X = NspksAll(:, [nn nnn]);
%                         Y = FF;
%                         
%                         if tempfix ==1
%                             if size(X,1)<5
%                                 continue
%                             end
%                         end
%                         
%                         
%                         if handcodeinteraction==1
%                             % ---- make new column by multiplying X (and then subtrctig mean)
%                             % ---- since other way (using interactions) can
%                             % have multiplication of negative numbers (if
%                             % subtract mean first).
%                             X = [X, X(:,1).* X(:,2)];
%                             modelspec = 'linear';
%                         else
%                             modelspec = 'interactions';
%                         end
%                         
%                         % --- subtract mean
%                         X = X - repmat(mean(X,1), size(X,1), 1);
%                         Y = FF - mean(FF);
%                         
%                         
%                         mdl = fitlm(X, Y, modelspec);
%                         
%                         % ==== extract coeff, se, and p, from mdl
%                         betas = mdl.Coefficients.Estimate(2:end);
%                         pvals = mdl.Coefficients.pValue(2:end);
%                         betasSE = mdl.Coefficients.SE(2:end);
%                         
%                         % ========= collect
%                         All_NeurPairFF_beta = [All_NeurPairFF_beta; betas']; % trials x 3(2 regions + interactio)
%                         All_NeurPairFF_pval = [All_NeurPairFF_pval; pvals'];
%                         All_NeurPairFF_Brain = [All_NeurPairFF_Brain; bregtmp];
%                         All_NeurPairFF_NeurID = [All_NeurPairFF_NeurID; ...
%                             [neurons_thisset(nn) neurons_thisset(nnn)]];
%                         All_NeurPairFF_BirdID = [All_NeurPairFF_BirdID; ...
%                             i];
%                         All_NeurPairFF_ExptID = [All_NeurPairFF_ExptID; ...
%                             ii];
%                         All_NeurPairFF_Ntrials = [All_NeurPairFF_Ntrials; ...
%                             length(FF)];
%                         All_NeurPairFF_motifnum = [All_NeurPairFF_motifnum; ...
%                             mm];
%                         All_NeurPairFF_setnum = [All_NeurPairFF_setnum; ...
%                             iii];
%                         
%                         
%                         % ================= PLOT EXAMPLES
%                         if plotRaw==1
%                             % -- region 1 vs FF
%                             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                             title(motifstr);
%                             ylabel('FF');
%                             xlabel(['n' num2str(neurons_thisset(nn)) '-' BrainregionsList{nn}]);
%                             lt_regress(Y, X(:,1), 1);
%                             lt_plot_zeroline;
%                             lt_plot_zeroline_vert
%                             
%                             % -- region 2 vs. FF
%                             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                             title(motifstr);
%                             ylabel('FF');
%                             xlabel(['n' num2str(neurons_thisset(nnn)) '-' BrainregionsList{nnn}]);
%                             lt_regress(Y, X(:,2), 1);
%                             lt_plot_zeroline;
%                             lt_plot_zeroline_vert
%                             
%                             
%                             % -- region1*region2 vs. FF
%                             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                             title(motifstr);
%                             ylabel('FF');
%                             xlabel('interaction');
%                             lt_regress(Y, X(:,1).*X(:,2), 1);
%                             lt_plot_zeroline;
%                             lt_plot_zeroline_vert
%                             
%                             %                             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                             %                         title(motifstr);
%                             %                         lt_plot_stem3(X(:,1), X(:,2), Y, 'r', 1);
%                             %
%                             %                         xlabel(['n' num2str(neurons_thisset(nn)) '-' BrainregionsList{nn}]);
%                             %                         ylabel(['n' num2str(neurons_thisset(nnn)) '-' BrainregionsList{nnn}]);
%                             %                         zlabel('FF');
%                             %                         lt_plot_zeroline;
%                             %                             lt_plot_zeroline_vert;
%                         end
%                     end
%                 end
%                 
%                 if plotRaw==1
%                     pause
%                     close all;
%                 end
                
%                 %% 2) EACH NEURON REGRESS WITH FF
%                 
%                 neurFF_brain = BrainregionsList;
%                 
%                 if isnan(FF(1))
%                     All_NeurFF_Rho = [All_NeurFF_Rho nan(1,nneur)];
%                     All_NeurFF_P = [All_NeurFF_P nan(1,nneur)];
%                     All_NeurFF_Brain = [All_NeurFF_Brain neurFF_brain];
%                     
%                 else
%                     
%                     [rho, p] = corr([NspksAll FF]);
%                     neurFF_rho = rho(1:end-1,end); % corr between individual neurons and FF (in order)
%                     neurFF_p = p(1:end-1, end);
%                     
%                     All_NeurFF_Rho = [All_NeurFF_Rho neurFF_rho'];
%                     All_NeurFF_P = [All_NeurFF_P neurFF_p'];
%                     All_NeurFF_Brain = [All_NeurFF_Brain neurFF_brain];
%                     
%                 end
                
                %% 3) CORRELATION BETWEEN NEURONS
                for nn=1:nneur
                    for nnn=nn+1:nneur
                        
                        % ============================ SPIKE CORRELATIONS
%                         rtmp = rho(nn, nnn);
%                         ptmp = p(nn, nnn);
                        bregtmp = {BrainregionsList{nn}, BrainregionsList{nnn}};
                        
                        % ===== plot this as example?
                        if plotRaw_CC==1 & any(strcmp(bregtmp, 'LMAN')) & ...
                                any(strcmp(bregtmp, 'RA'))
                            plotThisCC = 1;
                        else
                            plotThisCC =0;
                        end
                        
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
                        
                        % ###################################### CALCULATE
                        ccRealAll = [];
                        ccShiftAll = [];
                        
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
                            
                            ccRealAll = [ccRealAll; cc'];
                           
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
                            
                            
                            % ################################# shift control
                            if t<ntrials
                                fr2_shift = dattmp2.SegmentsExtract(t+1).FRsmooth_rate_CommonTrialDur;
                            else
                                fr2_shift = dattmp2.SegmentsExtract(1).FRsmooth_rate_CommonTrialDur;
                            end
                            fr2shift_short = fr2_shift(frx>=windtmp(1) & frx<windtmp(2));
                            
                            
                            [cc, lags] = xcov(fr1_short, fr2shift_short, xcovwindmax/0.001, 'coeff');
                            
                            ccShiftAll = [ccShiftAll; cc'];
                            
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
                        
                        
                        % ============================ OUTPUT
                        MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).CCreal = ccRealAll;
                        MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).CCshift = ccShiftAll;
                        MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).x = 0.001*lags;
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
