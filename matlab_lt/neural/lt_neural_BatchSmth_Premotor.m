function lt_neural_BatchSmth_Premotor(DATAllSwitches, motifstoplot, premotor_wind, ...
    plotRaw, removeNoiseTrials)

%% lt 11/17/17 - plot premotor and change in premotor for all chans and blocks

% plotRaw = 0; % if 1, then plots all trials and blocks overlaied. if 0 then goes straight to cross corr summary.
% 
% % ====== input params
% motifstoplot = {'a(b)', 'j(b)', 'h(g)'}; % if empty, plots all
% premotor_wind = [-0.03 0.02]; % for cross correlation
% removeNoiseTrials = 1; if 1, then throws out all trials that are noise
% (based on [fname].noise.mat).


%%
RhoAllStruct = struct;
chanstoplot = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));
numswitches = length(DATAllSwitches.switch);
nummotifs = length(DATAllSwitches.switch(1).motif);


%%
for cc = chanstoplot
    for mm = 1:nummotifs
        
        motif = DATAllSwitches.switch(1).motif(mm).batchinorder(1).motifname;
        if ~isempty(motifstoplot)
            if ~any(ismember(motifstoplot, motif))
                continue
            end
        end
        
        figcount=1;
        subplotrows=4;
        subplotcols=6;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        DatMeansAllSwitches = []; % to plot timecourse summary
        for i=1:numswitches
            
            DatMeans ={};
            DatSEMs = {};
            CondAll = {};
            
            for bb = 1:length(DATAllSwitches.switch(i).motif(mm).batchinorder)
                % i.e., each batch corresponds to a "condition" (e.g.
                % WNon/off; or DIR/UNDIR);
                datmat = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAll{cc};
                t = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).t;
                FF = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).FF;
                cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                motif = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).motifname;
                fs = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).fs;
                pretime = DATAllSwitches.switch(i).params.pretime;
                pretime = pretime/fs;
                
                noisetrials = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).NoiseTrials{cc};
                noisetrials = logical(noisetrials);
                
                % ==================================
                if removeNoiseTrials ==1
                    datmat(noisetrials,:) = [];
                    FF(noisetrials) = [];
                end
                
                if isempty(datmat)
                    datmedian = nan(size(t));
                    datSEM =  nan(size(t));
                    
                else
                    % ======== plot all dats, and overlay median
                    if plotRaw==1
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title(['ch' num2str(cc) '-sw' num2str(i) '-' motif '[' cond ']']);
                        
                        plot(t, datmat, 'Color', [0.7 0.7 0.7]);
                    end
                    
                    % -- overlay median, with 75th tiles
                    datmedian = median(datmat,1);
                    if size(datmat,1)>1
                    datSEM = lt_sem(datmat);
                    else
                        datSEM = nan(size(t));
                    end
                    
                    if plotRaw==1 & size(datmat,1)>1
                        shadedErrorBar(t, datmedian, datSEM, {'Color','r'}, 1);
                        line([pretime pretime],ylim);
                        lt_plot_zeroline;
                        axis tight;
                    end
                 end   
                    % ============ save
                    DatMeans = [DatMeans datmedian];
                    DatSEMs = [DatSEMs datSEM];
                    CondAll = [CondAll cond];
                    
                    DatMeansAllSwitches = [DatMeansAllSwitches; datmedian];
                
            end
            
            % ================= PLOT EACH BATCH'S MEAN, OVERLAYED
            if plotRaw==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
            end
            
            for bb = 1:length(DatMeans)
                
                if strcmp(CondAll{bb}, 'WNon')
                    plotcol = 'r';
                elseif strcmp(CondAll{bb}, 'DIR')
                    plotcol = 'b';
                else
                    plotcol = 'k';
                end
                if plotRaw==1
                    shadedErrorBar(t, DatMeans{bb}, DatSEMs{bb}, {'Color', plotcol}, 1);
                    line([pretime pretime],ylim);
                    lt_plot_zeroline;
                    axis tight;
                end
                
            end
            
        end
        if plotRaw==1
            linkaxes(hsplots, 'xy');
        end
        
        % ===== plot timecourse showing how mean contour correlates with
        % that of the first batch
        % ---- 1) get premotor window
        assert(length(t) == size(DatMeansAllSwitches,2), 'asdfasd');
        windwind = premotor_wind+pretime;
        
        inds = t>=windwind(1) & ...
            t<=windwind(2);
        DatMeansAllSwitches = DatMeansAllSwitches(:, inds);
        t = t(inds);
        
        % ----- 2) bin to 1ms
        TrimDown = 1;
        binsize = 0.001;
        [DatMeansAllSwitches, t] = lt_neural_v2_QUICK_binFR(DatMeansAllSwitches, t, binsize, TrimDown);
        
        % ---------- each one get cross corr versus first on
        datfirst = DatMeansAllSwitches(1,:);
        RhoAll = [];
        for j=1:size(DatMeansAllSwitches,1)
            
            datthis = DatMeansAllSwitches(j,:);
            
            rho = corr(datfirst', datthis');
            RhoAll = [RhoAll rho];
        end
        if plotRaw==1
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel('rho, versus first block');
            xlabel('block num');
            plot(1:length(RhoAll), RhoAll, '-ok');
        end
        
        RhoAllStruct.chan(cc).motifs(mm).rhoacrossbatch = RhoAll;
        RhoAllStruct.chan(cc).motifs(mm).motif = motif;
        
        
        % ============================= EACH SWITCH GET CORR BETWEEN THE
        % PRE AND POST
        rho_prevspost_withinswitch = [];
        frdiff_postminuspre_withinswitch =[];
        
        for j=1:numswitches
           ind1 = j*2-1;
           ind2 = j*2;
           
           dat1 = DatMeansAllSwitches(ind1, :);
           dat2 = DatMeansAllSwitches(ind2, :);
           
           % ---- correlation, pre vs. post
           rhotmp = corr(dat1', dat2');
           rho_prevspost_withinswitch = [rho_prevspost_withinswitch rhotmp];
           
 
           % ---- mean difference, pre vs. post
           frdiff = mean(dat2) - mean(dat1);
           frdiff_postminuspre_withinswitch =[frdiff_postminuspre_withinswitch ...
               frdiff];
        end
        
        
        if plotRaw==1
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel('rho, pre vs. post');
            xlabel('switch num');
            plot(1:length(rho_prevspost_withinswitch), ...
                rho_prevspost_withinswitch, '-ok');
            xlim([0 length(rho_prevspost_withinswitch)+1]);
            ylim([-1 1]);
            lt_plot_zeroline;
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel('frdiff, post minus pre');
            xlabel('switch num');
            plot(1:length(frdiff_postminuspre_withinswitch), ...
                frdiff_postminuspre_withinswitch, '-ok');
            
            xlim([0 length(rho_prevspost_withinswitch)+1]);
            ylim([-20 20]);
            lt_plot_zeroline;
            
        end
        
        RhoAllStruct.chan(cc).motifs(mm).rho_prevspost_withinswitch = rho_prevspost_withinswitch;
        RhoAllStruct.chan(cc).motifs(mm).frdiff_postminuspre_withinswitch = frdiff_postminuspre_withinswitch;
        
    end
end


%% ============================== PLOT SUMMARY
% ------- each block (corr with baseline)
figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
motifsexist = find(~cellfun('isempty', {RhoAllStruct.chan(chanstoplot(1)).motifs.rhoacrossbatch}));
for i = 1:length(motifsexist)
    mm= motifsexist(i);
    
    plotcols = lt_make_plot_colors(length(chanstoplot), 0, 0);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([RhoAllStruct.chan(chanstoplot(1)).motifs(mm).motif]);
    xlabel('block num');
    ylabel('corr');
    for cc = 1:length(chanstoplot)
        chan = chanstoplot(cc);
        
        dat = RhoAllStruct.chan(chan).motifs(mm).rhoacrossbatch;
        plot(1:length(dat), dat, '-o', 'Color', plotcols{cc});
        lt_plot_text(length(dat), dat(end), num2str(chan), plotcols{cc});
    end
    
    ylim([-0.5 1]);
    lt_plot_zeroline;
    
end

% --------- each switch, pre vs. post corr;
figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots
motifsexist = find(~cellfun('isempty', {RhoAllStruct.chan(chanstoplot(1)).motifs.rhoacrossbatch}));
for i = 1:length(motifsexist)
    mm= motifsexist(i);
    
    plotcols = lt_make_plot_colors(length(chanstoplot), 0, 0);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([RhoAllStruct.chan(chanstoplot(1)).motifs(mm).motif])
    xlabel('switch num');
    ylabel('corr');
    hsplots = [hsplots hsplot];
    for cc = 1:length(chanstoplot)
        chan = chanstoplot(cc);
        
        dat = RhoAllStruct.chan(chan).motifs(mm).rho_prevspost_withinswitch;
        plot(1:length(dat), dat, '-o', 'Color', plotcols{cc});
        lt_plot_text(length(dat), dat(end), num2str(chan), plotcols{cc});
    end
    
    ylim([-0.5 1]);
    lt_plot_zeroline;
    
end
linkaxes(hsplots, 'xy');

% --------- each switch, post minus pre FR diff;
figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots =[];

motifsexist = find(~cellfun('isempty', {RhoAllStruct.chan(chanstoplot(1)).motifs.rhoacrossbatch}));
for i = 1:length(motifsexist)
    mm= motifsexist(i);
    
    plotcols = lt_make_plot_colors(length(chanstoplot), 0, 0);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([RhoAllStruct.chan(chanstoplot(1)).motifs(mm).motif])
    xlabel('switch num');
    ylabel('FR diff (post minus pre)');
    hsplots = [hsplots hsplot];
    for cc = 1:length(chanstoplot)
        chan = chanstoplot(cc);
        
        dat = RhoAllStruct.chan(chan).motifs(mm).frdiff_postminuspre_withinswitch;
        plot(1:length(dat), dat, '-o', 'Color', plotcols{cc});
        lt_plot_text(length(dat), dat(end), num2str(chan), plotcols{cc});
    end
    
    ylim([-20 20]);
    lt_plot_zeroline;
    
end
linkaxes(hsplots, 'xy');

