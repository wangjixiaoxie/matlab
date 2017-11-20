function lt_neural_BatchSmth_CtxtSep(DATAllSwitches, MotifSets, premotor_wind, ...
    useCorr, plotRaw)
%% lt 11/17/17 - smoothed FR separation between contexts for the same sylabel.


%% ====== input params
% clear MotifSets;
% MotifSets{1} = {'a(b)', 'j(b)'};
% MotifSets{2} = {'(a)ab', 'a(a)b'};
% MotifSets{3} = {'(j)jb', 'j(j)b'};
% MotifSets{4} = {'jb(h)', 'jbh(h)'};
% 
% plotRaw =0; % if 1, then plots FR traces (USEFUL). if 0, then just plots summary of corelation coeff.
% useCorr =1; % if 1, then cauclate pearson's corr, if 0, then euclid dist (5ms bins), to ask about similarity,
% 
% premotor_wind = [-0.03 0.02]; % for cross correlation

%%
nummotifpairs = length(MotifSets);
hsplots =[];
chanstoplot = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));
numswitches = length(DATAllSwitches.switch);
RhoAll_AcrossMotifSets = struct;


%%
for cc = chanstoplot
    
    RhoAll_AcrossMotifSets.chan(cc).motifset = {};
    for mm =1:length(MotifSets)
        figcount=1;
        subplotrows=4;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots1 = [];
        hsplots2 =[];
        
        RhoAll = []; % for each block, extracts average correlation coeff between premotor activity for motifs;
        for i=1:numswitches
            
            plotcols = lt_make_plot_colors(length(MotifSets{mm}), 0, 0);
            
            for bb = 1:length(DATAllSwitches.switch(i).motif(mm).batchinorder)
                
                if plotRaw==1
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots1 = [hsplots1 hsplot];
                    cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                    title(['ch' num2str(cc) '-sw' num2str(i) '[' cond ']']);
                end
                
                % ============= go thru all motifs in this motif set and
                % plot for this epoch
                motifstoplot = MotifSets{mm};
                FFvals = cell(1,length(motifstoplot));
                DatMedianAll = cell(1,length(motifstoplot));
                for nn = 1:length(motifstoplot)
                    
                    motif = motifstoplot{nn};
                    
                    ind = strcmp({DATAllSwitches.switch(i).motif.motif}, motif);
                    
                    assert(sum(ind)==1,'asdf');
                    
                    datmat = DATAllSwitches.switch(i).motif(ind).batchinorder(bb).DatAll{cc};
                    t = DATAllSwitches.switch(i).motif(ind).batchinorder(bb).t;
                    fs = DATAllSwitches.switch(i).motif(ind).batchinorder(bb).fs;
                    pretime = DATAllSwitches.switch(i).motif(ind).batchinorder(bb).pretime./fs;
                    % =========== plot median
                    datmedian = median(datmat,1);
                    datSEM = lt_sem(datmat);
                    
                    if plotRaw==1
                        shadedErrorBar(t, datmedian, datSEM, {'Color', plotcols{nn}}, 1);
                        lt_plot_text(t(50), datmedian(50)+5*nn, motif, plotcols{nn});
                        line([pretime pretime], ylim, 'Color','b');
                    end
                    
                    % ============= collect
                    FF = DATAllSwitches.switch(i).motif(ind).batchinorder(bb).FF;
                    FFvals{nn} = FF;
                    DatMedianAll{nn} = datmedian;
                end
                
                % ========================== PLOT FF VALS DISTRIBUTION
                if plotRaw==1
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots2 = [hsplots2 hsplot];
                    title(['ch' num2str(cc) '-sw' num2str(i) '[' cond ']']);
                    lt_plot_MultDist(FFvals, 1:length(FFvals), 1, 'k', 0, 0);
                    set(gca, 'XTickLabel', MotifSets{mm});
                    rotateXLabels(gca, 45);
                    axis tight;
                end
                
                % ============================= CORRELATION BETWEEN TRACES,
                %, ONE FOR EACH BATCH IN THE SWITCH
                % ------- only take premotor window;
                rhoall = [];
                for nn=1:length(DatMedianAll)
                    for pp = nn+1:length(DatMedianAll)
                        
                        windwind = premotor_wind+pretime;
                        
                        inds = t>=windwind(1) & ...
                            t<=windwind(2);
                        dtmp1 = DatMedianAll{nn}(inds);
                        dtmp2 = DatMedianAll{pp}(inds);
                        t = t(inds);
                        
                        % -------- bin
                            TrimDown = 1;
                            binsize = 0.001;
                            [dtmp1, xtimes] = lt_neural_v2_QUICK_binFR(dtmp1, t, binsize, TrimDown);
                            [dtmp2, xtimes] = lt_neural_v2_QUICK_binFR(dtmp2, t, binsize, TrimDown);
                            
                        % ===============================
                        if useCorr==1
                            rho = corr(dtmp1, dtmp2);
                            
                        elseif useCorr ==0
                            % then use euclidian distance                            
                            rho = sqrt(sum((dtmp2-dtmp1).^2));
                            
                        end
                        rhoall = [rhoall rho];
                    end
                end
                rhoall = mean(rhoall); % average over all pairs;
                RhoAll = [RhoAll rhoall];
            end
        end
        
        % =============== plot average correlation over blocks
        if plotRaw==1
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel('mean corr between motifs, single block');
            xlabel('block num');
            lt_plot_zeroline;
            ylim([-0.2 1]);
            plot(1:length(RhoAll), RhoAll, '-ok');
            
            linkaxes(hsplots1, 'xy');
            linkaxes(hsplots2, 'xy');
        end
        RhoAll_AcrossMotifSets.chan(cc).motifset = [RhoAll_AcrossMotifSets.chan(cc).motifset ...
            RhoAll];
        
                
    end
end


% =============================================== PLOT ALL CORRELATION
% COEFFICIENTS
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
for mm=1:length(MotifSets)
    
    plotcols = lt_make_plot_colors(length(chanstoplot), 0, 0);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([MotifSets{mm}])
    for cc = 1:length(chanstoplot)
        chan = chanstoplot(cc);
        
        dat = RhoAll_AcrossMotifSets.chan(chan).motifset{mm};
        plot(1:length(dat), dat, '-o', 'Color', plotcols{cc});
        lt_plot_text(length(dat), dat(end), num2str(chan), plotcols{cc});
    end
    if useCorr==1
    ylim([-0.5 1]);
    end
    lt_plot_zeroline;
    
end