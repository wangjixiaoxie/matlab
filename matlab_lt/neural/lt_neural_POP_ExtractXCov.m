function MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk)
%% extracts pairwise cross correlations. sticks back into structure

%% PARAMS
onlySetsWithRALMANpaired = 0;
plotRaw_CC = 0;
plotSummary = 0;

plotCols = lt_make_plot_colors(3, 0, 0);
plotColsRegions = {'LMAN', 'RA', 'X'};

warnedflag = 0;

if ~exist('xcov_dattotake', 'var')
    xcov_dattotake = []; % leave empty to use all data in
    % xcov_dattotake = [-0.1 0.05]; % rel syl onset; % leave empty to take all
end

% ========== xcorr analysis of frate
if isempty(xcovwindmax)
    xcovwindmax = 0.05; %
end

% ====== binning spikes:
if ~exist('binsize_spk', 'var')
    binsize_spk = 0.005; % default, 5ms bins for cross corr
end


%% RUNS

NumBirds = length(MOTIFSTATS_pop.birds);

for i=1:NumBirds
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
    
    for ii=1:numexpts
        exptname = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
        for iii=1:numsets
            disp(['bird' num2str(i) '-expt' num2str(ii) '-set' num2str(iii)])
            
            neurons_thisset = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{iii};
            nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif);
            nneur = length(neurons_thisset);
            
            % =========== WHAT BRAIN REGIONS ARE THESE?
            BrainregionsList = {SummaryStruct.birds(i).neurons(neurons_thisset).NOTE_Location};
            
            
            
            %% ######################### LMAN/RA PAIRS [MODIFY TO BE GENERAL]
            if onlySetsWithRALMANpaired==1
                
                indsLMAN = strcmp(BrainregionsList, 'LMAN');
                indsRA = strcmp(BrainregionsList, 'RA');
                indsX = strcmp(BrainregionsList, 'X');
                
                if ~any(indsLMAN) & ~any(indsRA)
                    continue
                end
            end
            
            %% GO THRU MOTIFS AND COLLECT THINGS
            
            for mm = 1:nummotifs
                
                segextract_for_trialdur = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                motifstr = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).regexpstr;
                motifpredur = MOTIFSTATS_pop.birds(i).params.motif_predur;
                motifpostdur =  MOTIFSTATS_pop.birds(i).params.motif_postdur;
                % ========== for plotting raw
                figcount=1;
                subplotrows=5;
                subplotcols=3;
                fignums_alreadyused=[];
                hfigs=[];
                
                
                
                %% =================== CORRELATION BETWEEN NEURONS
                for nn=1:nneur
                    for nnn=nn+1:nneur
                        
                        
                        bregtmp = {BrainregionsList{nn}, BrainregionsList{nnn}};
                        MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).bregtmp = bregtmp;

                        % ===== plot this as example?
                        if plotRaw_CC==1 & any(strcmp(bregtmp, 'LMAN')) & ...
                                any(strcmp(bregtmp, 'RA'))
                            plotThisCC = 1;
                        else
                            plotThisCC =0;
                        end
                        
                        % ############################## CROSS CORRELATIONS
                        % ============================ EXTRACT DATA
                        dattmp1 = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(nn);
                        dattmp2 = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(nnn);
                        
                        clustnum1 = SummaryStruct.birds(i).neurons(neurons_thisset(nn)).clustnum;
                        clustnum2 = SummaryStruct.birds(i).neurons(neurons_thisset(nnn)).clustnum;
                        
                        
                        % =========================== EXTRACT FR AND SPIKES
                        % ----------- 1) SMOOTHED FR
                        dattmp1.SegmentsExtract = lt_neural_SmoothFR(dattmp1.SegmentsExtract, ...
                            clustnum1, '', '', '', segextract_for_trialdur);
                        dattmp2.SegmentsExtract = lt_neural_SmoothFR(dattmp2.SegmentsExtract, ...
                            clustnum2, '', '', '', segextract_for_trialdur);
                        
                        % --------------- 2) BIN SPIKES
                        % what is maximum common trial dur?
                        if (0)
                            stimes = {segextract_for_trialdur.spk_Times};
                        maxtime = cellfun(@max, stimes);
                        maxtime = max(maxtime);
                        else
                            if warnedflag==0
                                warnedflag=1
                            disp('NOTE: assuming is linear time warped, and max duration is relative to offset of motif');
                            pause
                            end
                           maxtime = segextract_for_trialdur(1).motifsylOnsets(end)+motifpostdur;
                        end
                        
                        assert(length(unique(cellfun(@unique, {dattmp1.SegmentsExtract.spk_Clust})))==1, 'other clusts here ...');
                        assert(length(unique(cellfun(@unique, {dattmp2.SegmentsExtract.spk_Clust})))==1, 'other clusts here ...');
                        
                        % Extract
                        [dattmp1.SegmentsExtract] = lt_neural_QUICK_SpkBinned(dattmp1.SegmentsExtract, maxtime, binsize_spk);
                        [dattmp2.SegmentsExtract] = lt_neural_QUICK_SpkBinned(dattmp2.SegmentsExtract, maxtime, binsize_spk);
                        
                     %% ================== CROSS-CORRELATION MATRIX BETWEEN NEURONS PAIRS
%                         if strcmp(exptname, 'RALMANlearn1') & nn==1 ...
%                                 & nnn==3
%                             keyboard
%                         end
%                             
                        
                        spkbin1 = [dattmp1.SegmentsExtract.spk_Binned];
                        spkbin2 = [dattmp2.SegmentsExtract.spk_Binned];
                        spkbin1 = spkbin1';
                        spkbin2 = spkbin2';
                        x = dattmp1.SegmentsExtract(1).spk_Binned_x;
                        
                        % -- auto1
                        rho1 = corr(double(spkbin1), double(spkbin1));
                                                
                        % -- auto2
                        rho2 = corr(double(spkbin2), double(spkbin2));
                        
                        % -- cross
                        rho12 = corr(double(spkbin1), double(spkbin2));
                        
                        if (0)
                           lt_figure; hold on;
                        
                           % --- cross
                           lt_subplot(2,2,1); hold on;
                           title('cross');
                           xlabel(bregtmp{2})
                           ylabel(bregtmp{1});
                           imagesc(x, x, rho12);
                           axis tight;
                           % --- put syl onsets offsets
                           ontimes = segextract_for_trialdur.sylOnTimes_RelDataOnset;
                           offtimes = segextract_for_trialdur.sylOffTimes_RelDataOnset;
                           
                           for lll =1:length(ontimes)

                               on = ontimes(lll);
                               off = offtimes(lll);
                               
                               line([on on], ylim, 'Color', 'g');
                               line([off off], ylim, 'Color', 'r');
                               line(xlim, [on on], 'Color', 'g');
                               line(xlim, [off off],'Color', 'r');
                           end
                           
                            % --- auto1
                           lt_subplot(2,2,2); hold on;
                           title('auto');
                           xlabel(bregtmp{1})
                           ylabel(bregtmp{1});
                           imagesc(x, x, rho1);
                           axis tight;
                           % --- put syl onsets offsets
                           ontimes = segextract_for_trialdur.sylOnTimes_RelDataOnset;
                           offtimes = segextract_for_trialdur.sylOffTimes_RelDataOnset;
                           
                           for lll =1:length(ontimes)

                               on = ontimes(lll);
                               off = offtimes(lll);
                               
                               line([on on], ylim, 'Color', 'g');
                               line([off off], ylim, 'Color', 'r');
                               line(xlim, [on on], 'Color', 'g');
                               line(xlim, [off off],'Color', 'r');
                           end
                          
                             % --- auto2
                           lt_subplot(2,2,3); hold on;
                           title('auto');
                           xlabel(bregtmp{2})
                           ylabel(bregtmp{2});
                           imagesc(x, x, rho2);
                           axis tight;
                           % --- put syl onsets offsets
                           ontimes = segextract_for_trialdur.sylOnTimes_RelDataOnset;
                           offtimes = segextract_for_trialdur.sylOffTimes_RelDataOnset;
                           
                           for lll =1:length(ontimes)

                               on = ontimes(lll);
                               off = offtimes(lll);
                               
                               line([on on], ylim, 'Color', 'g');
                               line([off off], ylim, 'Color', 'r');
                               line(xlim, [on on], 'Color', 'g');
                               line(xlim, [off off],'Color', 'r');
                           end
                          
                        end
                        
                        
                      %% ================== CROSS CORRELATION [uses one of 2 methods (fr or spkbins)
                        % -------------------- PARAMS, INITIATE
                        ntrials = length(dattmp1.SegmentsExtract);
                        
                        ccRealAll = [];
                        ccPSTH = [];
                        ccShiftAll =[];
                        
                        ccAuto1 = [];
                        ccAuto1Shift = [];
                        
                        ccAuto2 = [];
                        ccAuto2Shift = [];
                        
                        if plotThisCC==1
                            lt_figure; hold on;
                            ypos = 1;
                            
                            lt_subplot(8, 2, 1:14); hold on;
                            title([bregtmp{1} '-' bregtmp{2} '[' motifstr ']-set' ...
                                num2str(iii) '-n' num2str(neurons_thisset(nn)) ',' ...
                                num2str(neurons_thisset(nnn))]);
                        end
                        
                        
                        if (0)
                            % ############## OLD VERSION - cross covariaince, uses shift
                            % predictor, and smoothe FR [CURRENTLY NOT
                            % OUTPUTING...]
                            for t = 1:ntrials
                                
                                fr1 = dattmp1.SegmentsExtract(t).FRsmooth_rate_CommonTrialDur;
                                fr2 = dattmp2.SegmentsExtract(t).FRsmooth_rate_CommonTrialDur;
                                frx = dattmp1.SegmentsExtract(t).FRsmooth_xbin_CommonTrialDur;
                                
                                % --- limit to data in "premotor" window
                                if isempty(xcov_dattotake)
                                    windtmp = [-100000 1000000]; % crrazy large
                                else
                                    windtmp = motifpredur + xcov_dattotake;
                                end
                                
                                fr1_short = fr1(frx>=windtmp(1) & frx<windtmp(2));
                                fr2_short = fr2(frx>=windtmp(1) & frx<windtmp(2));
                                
                                % ------------------- calculate
                                [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001);
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
                                    fr2_shift = dattmp2.SegmentsExtract(t-1).FRsmooth_rate_CommonTrialDur;
                                end
                                
                                fr2shift_short = fr2_shift(frx>=windtmp(1) & frx<windtmp(2));
                                
                                
                                [cc, lags] = xcov(fr1_short, fr2shift_short, xcovwindmax/0.001');
                                
                                ccShiftAll = [ccShiftAll; cc'];
                            end
                            
                            % ---- just to plot
                            %                             ccPSTH = mean(ccShiftAll,1)';
                            binsize_spk = 0.001;
                        else
                            %% ============== NEW VERSION - bins spikes, cross-correlation,
                            % corrected against cross corr of PSTH
                            for t = 1:ntrials
                                
                                % -------- extract spks binned
                                spkbin1 = dattmp1.SegmentsExtract(t).spk_Binned;
                                spkbin2 = dattmp2.SegmentsExtract(t).spk_Binned;
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp1.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin1 = spkbin1(indtmp);
                                    spkbin2 = spkbin2(indtmp);
                                end
                                
                                
                                % ------------------- calculate
                                %                                 [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
                                [cc, lags] = xcorr(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                
                                ccRealAll = [ccRealAll; cc'];
                                
                                
                                % ========================== SHIFT
                                % PREDICTOR
                                if t<ntrials
                                    spkbin2shift = dattmp2.SegmentsExtract(t+1).spk_Binned;
                                else
                                    spkbin2shift = dattmp2.SegmentsExtract(t-1).spk_Binned;
                                end
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp1.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin2shift = spkbin2shift(indtmp);
                                end
                                
                                % ------------------------------ calculate
                                [cc, lags] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                
                                ccShiftAll = [ccShiftAll; cc'];
                            end
                            
                            %% ========================= AUTO COVARIANCE (1)
                            for t = 1:ntrials
                                
                                % -------- extract spks binned
                                spkbin1 = dattmp1.SegmentsExtract(t).spk_Binned;
                                spkbin2 = dattmp1.SegmentsExtract(t).spk_Binned;
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp1.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin1 = spkbin1(indtmp);
                                    spkbin2 = spkbin2(indtmp);
                                end
                                
                                
                                % ------------------- calculate
                                %                                 [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
                                [cc, lags] = xcorr(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                
                                ccAuto1 = [ccAuto1; cc'];
                                
                                
                                % ========================== SHIFT
                                % PREDICTOR
                                if t<ntrials
                                    spkbin2shift = dattmp1.SegmentsExtract(t+1).spk_Binned;
                                else
                                    spkbin2shift = dattmp1.SegmentsExtract(t-1).spk_Binned;
                                end
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp1.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin2shift = spkbin2shift(indtmp);
                                end
                                
                                % ------------------------------ calculate
                                [cc, lags] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                
                                ccAuto1Shift = [ccAuto1Shift; cc'];
                            end
                            
                         %% ========================= AUTO COVARIANCE (2)
                            for t = 1:ntrials
                                
                                % -------- extract spks binned
                                spkbin1 = dattmp2.SegmentsExtract(t).spk_Binned;
                                spkbin2 = dattmp2.SegmentsExtract(t).spk_Binned;
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp2.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin1 = spkbin1(indtmp);
                                    spkbin2 = spkbin2(indtmp);
                                end
                                
                                
                                % ------------------- calculate
                                %                                 [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
                                [cc, lags] = xcorr(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2, xcovwindmax/binsize_spk);
                                
                                ccAuto2 = [ccAuto2; cc'];
                                
                                
                                % ========================== SHIFT
                                % PREDICTOR
                                if t<ntrials
                                    spkbin2shift = dattmp2.SegmentsExtract(t+1).spk_Binned;
                                else
                                    spkbin2shift = dattmp2.SegmentsExtract(t-1).spk_Binned;
                                end
                                
                                % --- limit to data in "premotor" window
                                if ~isempty(xcov_dattotake)
                                    % shorten to desired window
                                    windtmp = motifpredur + xcov_dattotake;
                                    x = dattmp2.SegmentsExtract(t).spk_Binned_x;
                                    
                                    indtmp = x>=windtmp(1) & x<windtmp(2);
                                    spkbin2shift = spkbin2shift(indtmp);
                                end
                                
                                % ------------------------------ calculate
                                [cc, lags] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
                                
                                ccAuto2Shift = [ccAuto2Shift; cc'];
                            end

%%
                            
                            % ========================== pSTH control
                            spkall1 = [dattmp1.SegmentsExtract.spk_Binned];
                            spkall2 = [dattmp2.SegmentsExtract.spk_Binned];
                            
                            psth1 = mean(spkall1,2);
                            psth2 = mean(spkall2,2);
                            
                            
                            % --- limit to data in "premotor" window
                            if ~isempty(xcov_dattotake)
                                % shorten to desired window
                                windtmp = motifpredur + xcov_dattotake;
                                x = dattmp1.SegmentsExtract(t).spk_Binned_x;
                                indtmp = x>=windtmp(1) & x<windtmp(2);
                                psth1 = psth1(indtmp);
                                psth2 = psth2(indtmp);
                            end
                            
                            [ccPSTH, lags] = xcorr(psth1, psth2, xcovwindmax/binsize_spk);
                            %                             [ccPSTH, lags] = xcov(psth1, psth2, xcovwindmax/binsize_spk);
                            
                            % ################################## sanity
                            % check
                            if ~isempty(xcov_dattotake)
                                % then should all be same duration
                                assert(length(spkbin1) == length(spkbin2shift), 'sadf');
                                assert(length(spkbin2) == length(psth1), 'asfasd');
                            end
                            
                            
                            % ################################## OUTPUT
                            lags_sec = lags*binsize_spk;
                            
                            %% OUTPUT
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccAuto2 = ccAuto2;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccAuto2Shift = ccAuto2Shift;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccAuto1 = ccAuto1;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccAuto1Shift = ccAuto1Shift;
                            
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccRealAll = ccRealAll;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccShiftAll = ccShiftAll;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).ccPSTH = ccPSTH;
                            MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn).x = lags_sec;
                            
                            
                            %% ================ sanity check - plot CCs
                            if plotSummary==1 & any(strcmp(bregtmp, 'LMAN')) & any(strcmp(bregtmp, 'RA'))
                                lt_figure; hold on;
                                
                                % --------- 1) firing rates
                                lt_subplot(2,1,1); hold on;
                                title(motifstr);
                                ylabel([bregtmp{1} '(' num2str(neurons_thisset(nn)) ...
                                    ')-' bregtmp{2} '(' num2str(neurons_thisset(nnn)) ')']);
                                
                                y = [dattmp1.SegmentsExtract.FRsmooth_rate_CommonTrialDur];
                                x = 0.001*(1:size(y,1));
                                plot(x, y, '-', 'Color', [0.7 0.7 0.7]);
                                plot(x, mean(y'), 'k', 'LineWidth', 2);
                                ymax = max(y(:));
                                
                                y = [dattmp2.SegmentsExtract.FRsmooth_rate_CommonTrialDur];
                                x = 0.001*(1:size(y,1));
                                plot(x, y+1.2*ymax, '-', 'Color', [0.8 0.3 0.3]);
                                plot(x, mean(y')+1.2*ymax, 'r', 'LineWidth', 2);
                                
                                axis tight
                                
                                % --- lines for syl onsets, offsets
                                for kkk = 1:length(segextract_for_trialdur(1).motifsylOnsets)
                                    %                                     line([segextract_for_trialdur(1).motifsylOnsets(kkk) ...
                                    %                                         segextract_for_trialdur(1).motifsylOnsets(kkk)], ylim, 'Color', 'b');
                                    %                                     line([segextract_for_trialdur(1).motifsylOffsets(kkk) ...
                                    %                                         segextract_for_trialdur(1).motifsylOffsets(kkk)], ylim, 'Color', 'r');
                                    line([segextract_for_trialdur(1).motifsylOnsets(kkk) ...
                                        segextract_for_trialdur(1).motifsylOffsets(kkk)], [0 0], 'Color', 'b', ...
                                        'LineWidth', 2);
                                end
                                
                                
                                lt_subplot(2,2,3); hold on;
                                title('r,dash=PSTH; b,dash=SHIFT');
                                plot(lags*binsize_spk, ccRealAll', 'k-');
                                plot(lags*binsize_spk, mean(ccRealAll,1), 'r', 'LineWidth', 2);
                                plot(lags*binsize_spk, ccPSTH, 'r--', 'LineWidth', 2);
                                plot(lags*binsize_spk, mean(ccShiftAll,1), 'b--', 'LineWidth', 2);
                                axis tight;
                                lt_plot_zeroline;
                                lt_plot_zeroline_vert;
                                
                                lt_subplot(2,2,4); hold on;
                                plot(lags*binsize_spk, mean(ccRealAll,1) - ccPSTH', 'r-');
                                plot(lags*binsize_spk, mean(ccRealAll,1) - mean(ccShiftAll,1), 'b-');
                                axis tight;
                                lt_plot_zeroline;
                                lt_plot_zeroline_vert;
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
                        
                            
                        end
                    end
                    
                    
                end
            end
        end
    end
end
