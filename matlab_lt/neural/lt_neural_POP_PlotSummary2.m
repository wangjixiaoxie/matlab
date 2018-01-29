function lt_neural_POP_PlotSummary2(MOTIFSTATS_pop, SummaryStruct, DATSTRUCT)
%% lt 1/22/18 - plots output of lt_neural_POP_PlotSummary

%% PARAMS

smoothBeforeSubtr = 0; % smooth before subtract psth or shuffle corrector - see Kass


%% PLOT EACH MOTIF

bregionWanted = 'LMAN-RA';

numbirds = max(DATSTRUCT.Pairs.birdID);
numexpts = max(DATSTRUCT.Pairs.exptID);
nummotifs = max(DATSTRUCT.Pairs.motifnum);

figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


Ccallall = [];
for i=1:numbirds
    for ii=1:numexpts
        for mm =1:nummotifs
            
            inds = find(DATSTRUCT.Pairs.birdID==i & DATSTRUCT.Pairs.exptID==ii ...
                & DATSTRUCT.Pairs.motifnum==mm ...
                & strcmp(DATSTRUCT.Pairs.bregion_string	, bregionWanted));
            
            if isempty(inds)
                continue
            end
            
            
            % #################### PLOT EACH PAIR SEPARATELY
            CCfinalAll = [];
            for jj = inds'
                motifstr = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(1).motif(mm).regexpstr;
                bname = MOTIFSTATS_pop.birds(i).birdname;
                exptname = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
                neurons = DATSTRUCT.Pairs.neurIDreal(jj,:);
                
                % ------------------ 0) PLOT CROSS-COV
                ccRealAll = DATSTRUCT.Pairs.ccRealAll{jj};
                ccShiftAll = DATSTRUCT.Pairs.ccShiftAll{jj};
                x = DATSTRUCT.Pairs.xlags(jj,:);
                
                % --- get means
                ccreal = mean(ccRealAll,1);
                ccreal_sem = lt_sem(ccRealAll);
                ccneg = mean(ccShiftAll,1);
                ccneg_sem = lt_sem(ccShiftAll);
                
                ccreal_minus =  ccreal - ccneg;
                %
                %                 ccreal_minus_all = ccRealAll - repmat(ccneg, size(ccRealAll,1),1);
                %
                %% ==========
                if smoothBeforeSubtr==1
                    % -------- get gaussian window
                    filtdur = 0.015;
                    binsize = x(2)-x(1);
                    N = ceil(filtdur/binsize);
                    if mod(N,2)==0
                        N = N+1;
                    end
                    wind = gausswin(N);
                    wind = wind./sum(wind);
                    
                    % ------- smooth
                    ccreal_sm = conv(ccreal, wind);
                    edge = (length(ccreal_sm)-length(ccreal))/2;
                    ccreal_sm = ccreal_sm(edge+1:end-edge);
                    
                    ccneg_sm = conv(ccneg, wind);
                    edge = (length(ccneg_sm)-length(ccneg))/2;
                    ccneg_sm = ccneg_sm(edge+1:end-edge);
                    
                    % ------- get new final
                    ccreal_minus = ccreal_sm - ccneg_sm;
                    
                    if (0)
                        lt_figure; hold on ;
                        plot(x, ccreal, 'k');
                        plot(x, ccreal_sm, 'r');
                    end
                end
                
                %% ======== auto covariance
                if (0)
                ccAuto1minus = DATSTRUCT.Pairs.ccAuto1_minus(jj, :);
                ccAuto2minus = DATSTRUCT.Pairs.ccAuto2_minus(jj, :);
                
                ccAuto1 = DATSTRUCT.Pairs.ccAuto1_mean(jj,:);
                ccAuto2 = DATSTRUCT.Pairs.ccAuto2_mean(jj,:);
               
                
                ccreal_minus = ccreal./sqrt(ccAuto1.*ccAuto2);
                end
                %%
                
                % ---- plot mean raw vs. shuffle
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([bname '-' exptname '-' motifstr '-nn' num2str(neurons)]);
                shadedErrorBar(x, ccneg, ccneg_sem, {'Color', 'r', 'LineStyle', '--'}, 1)
                shadedErrorBar(x, ccreal, ccreal_sem, {'Color', 'k'}, 1)
                axis tight;
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                % ---- plot autocorr
                try
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title('auto, 1(k) and 2(b)');
                plot(x, ccAuto1, 'k');
                plot(x, ccAuto2, 'b');
                catch err
                    
                end
                % ---- plot shuffle subtracted
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                %                 title(['b- ' ]);
                plot(x, ccreal_minus, 'b');
                axis tight;
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                
                
                % ------------------- 1) PLOT SMOOTHED FR
                
                % ------------------- 2) OVERLAY SYL SEGMENTS
                
                % ========================= COLLECT
                CCfinalAll = [CCfinalAll; ccreal_minus];
                Ccallall = [Ccallall; ccreal_minus];
                
            end
            
            % ================== PLOT ALL PAIRS FOR THIS MOTIF
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('all neuron pairs');
            
            plot(x, CCfinalAll, 'b');
            axis tight;
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
            
        end
    end
end

%% =========== plot all
lt_figure; hold on;
title('all');

plot(x, Ccallall, 'Color', [0.7 0.7 0.7]);
ymean = mean(Ccallall);
ysem = lt_sem(Ccallall);
shadedErrorBar(x, ymean, ysem, {'Color', 'r'},1);

axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;





