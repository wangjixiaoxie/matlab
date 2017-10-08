function lt_neural_v2_CTXT_BranchEachPlot(ALLBRANCH, birdtoplot)

%%

motifpredur = ALLBRANCH.alignpos(1).ParamsFirstIter.motifpredur;
motifpostdur = ALLBRANCH.alignpos(1).ParamsFirstIter.motifpostdur;

%% First get dprime for all branches

% DPRIME STUFF
Niter = 3; % shuffles
Nmin = 3; % min sample size;
% Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;

% DprimeNegVersion = 'Wohl'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.

ALLBRANCH = lt_neural_v2_CTXT_AllBranchDprime(ALLBRANCH, Nmin, Niter, ...
    DprimeNegVersion);


%% ==== go thru each bird, plot all neurons and branches

numalignpos = length(ALLBRANCH.alignpos);

for i=1:numalignpos
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    
    
    for ii=1:numbirds
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        
        numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        figcount=1;
        subplotrows=4;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        if ~isempty(birdtoplot)
            if ~strcmp(birdname, birdtoplot)
                continue
            end
        end
        
        for iii=1:numbranches
            
            numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron);
            
            for nn=1:numneurons
                
                hsplots = [];
                sylregexp = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).prms_regexpstr;
                syllist = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).prms_regexpstrlist;
                dat = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn);
                
                if isempty(dat.xtimes)
                    continue
                end
                
                % ################################ fig 1 - raw hz
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' sylregexp '-n' num2str(nn)]);
                ylabel('fr(hz) OR classperformance/dprime*100');
                
                numclasses = length(dat.FR.classnum);
                plotcols = lt_make_plot_colors(numclasses, 0, 0);
                
                % ============ 1) mean FR
                for cc =1 :numclasses
                    
                    frmat = [dat.FR.classnum(cc).FRsmooth_rate_CommonTrialDur];
                    frmean = mean(frmat,2);
                    frsem = lt_sem(frmat');
                    x = dat.FR.classnum(cc).FRsmooth_xbin_CommonTrialDur-motifpredur;
                    shadedErrorBar(x, frmean, frsem, {'Color', plotcols{cc}}, 1);
                    
                    lt_plot_text(x(end), frmean(end), syllist{cc}, plotcols{cc});
                end
                
                
                % ============ 2) Syl contours
                for cc =1 :numclasses
                    
                    sylmean = dat.SylContoursByClass_means(cc,:);
                    sylstd = dat.SylContoursByClass_std(cc,:);
                    x = (1:length(sylmean))./1000;
                    x = x - motifpredur;
                    
                    shadedErrorBar(x, 20*sylmean-20, 20*sylstd, {'Color', plotcols{cc}}, 1);
                    
                end
                
                
                % ---------------------------
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                ylim([-30 200]);
                xlim([-motifpredur motifpostdur]);
                
                % ############### fig 2 - classifier and dprime
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' sylregexp '-n' num2str(nn)]);
                ylabel('classperformance (thick), dprime (thin)');
                
                % ============ 3) dprime
                % -- dat
                dprimemean = mean(100*dat.DprimeAllPairwise, 2);
                x = (1:length(dprimemean))./1000;
                x = x-motifpredur;
                
                plot(x, dprimemean, '-k');
                
                % -- neg
                dprimemean = mean(100*dat.DprimeAllPairwise_Neg, 2);
                x = (1:length(dprimemean))./1000;
                x = x-motifpredur;
                
                plot(x, dprimemean, '-r');
                
                
                % ============ 4) class performance
                lt_plot(dat.xtimes, 100*dat.yvals, {'LineStyle', '-', 'Color', 'k'});
                lt_plot(dat.xtimes, 100*dat.yvals_neg,{'LineStyle', '-', 'Color', 'r'});
                
                
                % ---------------------------
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                axis tight;
                xlim([-motifpredur motifpostdur]);
                
                linkaxes(hsplots, 'x');
                
                
            end
            
            
        end
        
        if isempty(birdtoplot)
        disp('PAUSSED --');
        pause;
        close all;
        end
    end
    
end