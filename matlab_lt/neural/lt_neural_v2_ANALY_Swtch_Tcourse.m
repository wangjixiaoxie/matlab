function lt_neural_v2_ANALY_Swtch_Tcourse(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get)
%% TO DO:
% LEARN FOR FOR NONNTARGET ARE STILL NOT IN CORRECT DIR - HAVE TO TAKE
% INTOA ACCOUNT MULTIPLE TARGETS

%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);

WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

FFsmthbinsize = 10;

minrends = 5; % for both train and base

premotorWind = SwitchStruct.params.premotorWind;


%% 

                    
                    if mod(FFsmthbinsize,2)==0
                        FFsmthbinsize= FFsmthbinsize+1; % conver to odd, so that median tval is at an actual datapoint.
                    end

%%
for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    if ~isempty(birdname_get)
        if ~strcmp(birdname, birdname_get)
            continue
        end
    end
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        if ~isempty(exptname_get)
            if ~strcmp(exptname, exptname_get)
                continue
            end
        end
        
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        targsyls = MotifStats.params.TargSyls;
        nummotifs = length(motiflist);
        
        WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
            MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            if ~isempty(switchnum_get)
                if switchnum_get ~= iii;
                    continue
                end
            end
            
            figcount=1;
            subplotrows=5;
            subplotcols=6;
            fignums_alreadyused=[];
            hfigs=[];
            
            
            goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] ...
                & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs]);
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
            swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            
            % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
            for j=1:nummotifs
                
                for nn=goodneurons
                    
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    
                    if ~isfield(segextract, 'FRsmooth_xbin')
                        continue
                    end
                    
                    baseInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds);
                    trainInds = find(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds);
                    
                    if length(baseInds)<minrends | length(trainInds) < minrends
                        % even for good neurosn, could occur if some motifs
                        % labeled pre but not post.
                        continue
                    end
                    
                    trialstoplot = [baseInds trainInds];
                    clustnum = MotifStats.neurons(nn).clustnum;
                    
                    
                    % ============= PLOT RASTER
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    
                    for tt = trialstoplot
                        indstmp = segextract(tt).spk_Clust==clustnum;
                        spktimes = segextract(tt).spk_Times(indstmp);
                        spktimes = spktimes(spktimes > WindowToPlot2(1) & ...
                            spktimes < WindowToPlot2(2));
                        for ttt =1:length(spktimes)
                            
                            line([spktimes(ttt) spktimes(ttt)], -[tt-0.4 tt+0.4], ...
                                'Color', plotcols{nn}, 'LineWidth', 1);
                        end
                    end
                    axis tight
                    set(gca, 'Ytick', []);
                    
                    % --- line for base vs. training
                    lastbaseind = max(baseInds);
                    line(xlim, [-lastbaseind-0.5 -lastbaseind-0.5], 'Color','k', 'LineWidth', 2);
                    
                    % -- line for syl onset
                    line([MotifStats.params.motif_predur MotifStats.params.motif_predur], ...
                        ylim, 'Color', 'k', 'LineWidth', 2);
                    
                    % --- line for neural analysis window
                    line([MotifStats.params.motif_predur+premotorWind(1) ...
                        MotifStats.params.motif_predur+premotorWind(1)], ylim,  ...
                        'LineStyle' ,'--', 'Color', 'k');
                    
                    line([MotifStats.params.motif_predur+premotorWind(2) ...
                        MotifStats.params.motif_predur+premotorWind(2)], ylim,  ...
                        'LineStyle' ,'--', 'Color', 'k');
                    
                    
                    % =============== PLOT SMOOTHED FR (BASE, TRAIN START,
                    % TRAIN END)
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    
                    % --- figure out max num rends
                    numtraintrials = length(trainInds);
                    maxnumrends = min([floor(numtraintrials/2), length(baseInds)]); % num base trials or 1/2 train trial, whichever is smaller
                    disp(maxnumrends);
                    
                    
                    % --- base
                    tbins = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    indstmp = max(baseInds)-maxnumrends+1:max(baseInds); % to equalize sample size
                    
                    smthFR = [segextract(indstmp).FRsmooth_rate_CommonTrialDur];
                    smthFR_mean = mean(smthFR, 2);
                    smthFR_sd = std(smthFR, 0, 2);
                    smthFR_sd = lt_sem(smthFR');
                    
                    indstmp = tbins>WindowToPlot2(1) & tbins<WindowToPlot2(2);
                    shadedErrorBar(tbins(indstmp), smthFR_mean(indstmp), smthFR_sd(indstmp), {'Color', [0.5 0.5 0.5]}, 1);
                    
                    % --- train start
                    indstmp = trainInds(1:maxnumrends);
                    
                    smthFR = [segextract(indstmp).FRsmooth_rate_CommonTrialDur];
                    smthFR_mean = mean(smthFR, 2);
                    smthFR_sd = std(smthFR, 0, 2);
                    smthFR_sd = lt_sem(smthFR');
                    
                    indstmp = tbins>WindowToPlot2(1) & tbins<WindowToPlot2(2);
                    shadedErrorBar(tbins(indstmp), smthFR_mean(indstmp), smthFR_sd(indstmp), {'Color', [0.7 0.3 0.3]}, 1);
                    
                    % --- train end
                    indstmp = trainInds(end-maxnumrends+1:end);
                    
                    smthFR = [segextract(indstmp).FRsmooth_rate_CommonTrialDur];
                    smthFR_mean = mean(smthFR, 2);
                    smthFR_sd = std(smthFR, 0, 2);
                    smthFR_sd = lt_sem(smthFR');
                    
                    indstmp = tbins>WindowToPlot2(1) & tbins<WindowToPlot2(2);
                    shadedErrorBar(tbins(indstmp), smthFR_mean(indstmp), smthFR_sd(indstmp), {'Color', [0.9 0.1 0.1]}, 1);
                    
                    
                    % -- line for syl onset
                    axis tight
                    line([MotifStats.params.motif_predur MotifStats.params.motif_predur], ...
                        ylim, 'Color', 'k', 'LineWidth', 2);
                    
                    % --- line for neural analysis window
                    line([MotifStats.params.motif_predur+premotorWind(1) ...
                        MotifStats.params.motif_predur+premotorWind(1)], ylim,  ...
                        'LineStyle' ,'--', 'Color', 'k');
                    
                    line([MotifStats.params.motif_predur+premotorWind(2) ...
                        MotifStats.params.motif_predur+premotorWind(2)], ylim,  ...
                        'LineStyle' ,'--', 'Color', 'k');
                    
                    % ================= PLOT TIMECOURSE OF NEURAL AND FF
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    
                    % ----- FF
                    ffvals = [segextract.FF_val];
                    tvals = [segextract.song_datenum];
                    % -- convert tvals to days
                    tmpday = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
                    tvals = lt_convert_EventTimes_to_RelTimes(datestr(tmpday, 'ddmmmyyyy'),...
                        tvals);
                    tvals = tvals.FinalValue;
                    tsmth = lt_running_stats(tvals, FFsmthbinsize);
                    
                    if ~all(isnan(ffvals));
                        % ---- get zscore
                        ffvals_basemean = mean(ffvals(baseInds));
                        ffvalsbaseSD = std(ffvals(baseInds));
                        
                        ffvals = (ffvals - ffvals_basemean)./ffvalsbaseSD;
                        
                        ffsmth = lt_running_stats(ffvals, FFsmthbinsize);
                        
                        plot(tsmth.Median, ffsmth.Median, 'kx');
                        %                       plot(tvals(trialstoplot), ffvals(trialstoplot), 'xk');
                        %
                    end
                    
                    % ---- neural (corr with baseline)
                    neuralmetricname = 'NEURvsbase_FRcorr';
                    neuralsim = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).(neuralmetricname);
                    %                    if ~strcmp(neuralmetricname, 'NEURvsbase_FRcorr')
                    % then get zscore
                    neurbasemean = mean(neuralsim(baseInds));
                    neurbaseSD = std(neuralsim(baseInds));
                    neuralsim = (neuralsim - neurbasemean)./neurbaseSD;
                    %                    end
                    neursmth = lt_running_stats(neuralsim, FFsmthbinsize);
                    %                     lt_plot(tsmth.Mean, neursmth.Mean, {'Errors', neursmth.SEM, ...
                    %                         'Color', plotcols{j}});
                    plot(tsmth.Median, neursmth.Mean, 'o', 'Color', plotcols{nn});
                    %                    plot(tvals(trialstoplot), neuralsim(trialstoplot), 'ob');
                    
                    
                    % --- stuff
                    axis tight;
                    ylim([-3 3]);
                    lt_plot_zeroline;
                    
                    % --- line for base vs. training
                    line([tvals(max(baseInds)) tvals(max(baseInds))], ylim, 'Color','k', 'LineWidth', 2);
                    
                    % --- title
                    if any(strcmp(targsyls, motiflist{j}))
                        title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}], 'Color', 'r');
                    else
                        title([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}]);
                    end
                end
            end
        end
    end
end




