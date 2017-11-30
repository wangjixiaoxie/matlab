function lt_neural_BatchSmth_CrossRegion(DATAllSwitches, motifstoplot,chanstoplot, windowmax, ...
    binsize, premotorwind, fs, plotRawDatOnly, removenoise)

%% 11/28/17 - LT

% motifstoplot = {'ab', 'jb', 'g'};
% motifstoplot = {'j(b)'};
% windowmax = 0.04;
% binsize = 0.002;
% premotorwind = [-0.035 0.025]; % use for ch14-21 on 11/12, morning.
% fs = 30000;
% plotRawDatOnly = 0;
% chanstoplot = [14 21];
% removenoise =1;

%%

nummotifs = length(DATAllSwitches.switch(1).motif);

%%
numswitches = length(DATAllSwitches.switch);

% --------- gaussian for smothing
if (0)
    windowsize=0.005; % from -2sd to +2sd
    sigma=(windowsize/4)*fs; %
    numsamps=4*sigma; % (get 2 std on each side)
    if mod(numsamps,2)==1
        numsamps = numsamps-1;
    end
    
    alpha= numsamps/(2*sigma); % N/2sigma
    gaussFilter = gausswin(numsamps, alpha);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
end

for k = 1:length(chanstoplot)
    chan1 = chanstoplot(k);
    
    for kk=k+1:length(chanstoplot)
        
        chan2 = chanstoplot(kk);
        
        for mm = 1:nummotifs
            
            motif = DATAllSwitches.switch(1).motif(mm).batchinorder(1).motifname;
            if ~isempty(motifstoplot)
                if ~any(ismember(motifstoplot, motif))
                    disp('SKIPPED');
                    continue
                end
            end
            
            % ====================== one figure for each motif
            figcount=1;
            if plotRawDatOnly==1
                subplotrows=5;
                subplotcols=3;
            else
                subplotrows=3;
                subplotcols=5;
            end
            fignums_alreadyused=[];
            hfigs=[];
            hsplots1 = [];
            hsplots2 = [];
            hsplots3 = [];
            
            for i=1:numswitches
                
                for bb = 1:length(DATAllSwitches.switch(i).motif(mm).batchinorder)
                    
                    % ================= for this channel pair, motif, and
                    % switch, get xcov in block preceding and following
                    % switch
                    cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                    %                     datmat1 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAllRaw{chan1};
                    %                     datmat2 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAllRaw{chan2};
                    datmat1 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAll{chan1};
                    datmat2 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAll{chan2};
                    t =  DATAllSwitches.switch(i).motif(mm).batchinorder(bb).t;
                    
                    
                    % ------------------ remove noise files
                    if removenoise ==1
                        noisetrials1 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).NoiseTrials{chan1};
                        noisetrials2 = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).NoiseTrials{chan2};
                        
                        % ---- union of noise trials
                        noisetrialsBoth = noisetrials1 | noisetrials2;
                        
                        datmat1(noisetrialsBoth, :) =[];
                        datmat2(noisetrialsBoth, :) = [];
                    end
                    
                    % ================= PLOT RAW DAT TO DETERMINE
                    % APPROPRIATE WINDOW TO TAKE TO AVOID NOISE
                    if plotRawDatOnly ==1
                        
                        numtrials = size(datmat1,1);
                        for tt=1:numtrials
                            
                            tTMP = (1/fs)*(1:size(datmat1,2));
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots3 = [hsplots3 hsplot];
                            title(['ch' num2str(chan1) '(k)vs' num2str(chan2) '(r), ' ...
                                motif ',sw' num2str(i) '[' cond ']']);
                            plot(tTMP, datmat1(tt,:), 'r');
                            plot(tTMP, datmat2(tt,:)+200, 'k');
                            axis tight
                        end
                        
                    else
                        
                        % ========== smooth rectify with small window
                        if (0)
                            datmat1 = abs(datmat1);
                            for nn =1:size(datmat1,1)
                                tmp = conv(datmat1(nn,:), gaussFilter);
                                tmp = tmp(:,numsamps/2:end-numsamps/2);
                                datmat1(nn,:) = tmp;
                                
                            end
                            %                         % -- clip off edges
                            %                         datmat1 = datmat1(:,numsamps/2:end-numsamps/2);
                            %
                            datmat2 = abs(datmat2);
                            for nn =1:size(datmat1,1)
                                tmp = conv(datmat2(nn,:), gaussFilter);
                                tmp = tmp(:,numsamps/2:end-numsamps/2);
                                datmat2(nn,:) = tmp;
                            end
                            
                        else
                            % this does nothing if using rectified
                            datmat1 = abs(datmat1);
                            datmat2 = abs(datmat2);
                        end
                        
                        if isempty(datmat1)
                            disp([motif '- SKIPPING! not enough tyrials'])
                            continue
                        end
                        
                        
                        % =================== cut off to premotor window
                        try
                            pretime = DATAllSwitches.switch(i).params.pretime;
                            pretime = pretime/fs;
                        catch err
                            disp('NOTE!!! making pretime 0.1');
                            pretime = 0.1;
                        end
                        
                        % ---------
                        tmp = pretime+premotorwind;
                        indstmp = t>=tmp(1) & t<=tmp(2);
                        datmat1 = datmat1(:, indstmp);
                        datmat2 = datmat2(:, indstmp);
                        t = t(indstmp);
                        assert(all(size(datmat1)==size(datmat2)), 'asdfasd');
                        
                        
                        % ---------------- bin activity
                        TrimDown = 1;
                        
                        [datmat1, xtimes] = lt_neural_v2_QUICK_binFR(datmat1, t, binsize, TrimDown);
                        [datmat2, xtimes] = lt_neural_v2_QUICK_binFR(datmat2, t, binsize, TrimDown);
                        
                        
                        % =================== calculate xcov
                        CCall = [];
                        numtrials = size(datmat2,1);
                        for tt=1:numtrials
                            [cc, lags] = xcov(datmat1(tt,:), datmat2(tt,:), ceil(windowmax/binsize), 'coeff');
                            CCall = [CCall; cc];
                        end
                        
                        % ================= calcualte xcov (shifted)
                        CCallSHUFFLE = [];
                        for tt=1:numtrials
                            
                            d1 = datmat1(tt,:);
                            if tt==numtrials
                                d2 = datmat2(1,:);
                            else
                                d2 = datmat2(tt+1,:);
                            end
                            
                            [cc, lags] = xcov(d1, d2, ceil(windowmax/binsize), 'coeff');
                            CCallSHUFFLE = [CCallSHUFFLE; cc];
                        end
                        
                        
                        % ===================== calculate noise correlation
                        frate1 = mean(datmat1,2);
                        frate2 = mean(datmat2,2);
                        windfrate = 10; % num trials
                        binfrate = 1; % trials shift
                        
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        xlabel(chan1); ylabel(chan2);
                        title('frate');
                        plot(frate1, frate2, 'ok');
                        axis tight;
                        
                        [ccfrate, lagfrate] = xcov(frate1, frate2, ceil(windfrate/binfrate), 'coeff');
                        
                        
                        
                        % ========================= PLOT ALL TRIALS [activity
                        % xcov]
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots1 = [hsplots1 hsplot];
                        title(['ch' num2str(chan1) 'vs' num2str(chan2) ', ' motif ',sw' num2str(i) '[' cond ']']);
                        plot(lags*binsize, CCall, 'Color', [0.7 0.7 0.7]);
                        lt_plot(lags*binsize, mean(CCall,1), {'Errors', lt_sem(CCall)});
                        lt_plot_zeroline;
                        
                        % --- overlay shifted
                        shadedErrorBar(lags*binsize, mean(CCallSHUFFLE,1), lt_sem(CCallSHUFFLE), ...
                            {'Color', [0.9 0.6 0.6]},1);
                        
                        
                        % ========================= PLOT ALL TRIALS [noise corr]
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots2 = [hsplots2 hsplot];
                        title('FRATE OVER TRIALS');
                        plot(lagfrate*binfrate, ccfrate, '-r');
                        ylim([-1 1]);
                        lt_plot_zeroline;
                    end
                end
            end
            if plotRawDatOnly==1
                linkaxes(hsplots3, 'xy');
            end
            linkaxes(hsplots1, 'xy');
            linkaxes(hsplots2, 'xy');
            pause;
            close all;
        end
        
    end
end

%%
if (0)
    % FOR PLOTTING DIAGNOSIS - comparing differeing smoothing and binning
    lt_figure; hold on;
    binsize = 0.001
    windowsize=0.005; % from -2sd to +2sd
    sigma=(windowsize/4)*fs; %
    numsamps=4*sigma; % (get 2 std on each side)
    if mod(numsamps,2)==1
        numsamps = numsamps-1;
    end
    
    alpha= numsamps/(2*sigma); % N/2sigma
    gaussFilter = gausswin(numsamps, alpha);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    
    DATMAT = datmat2;
    
    datmat1 = abs(DATMAT);
    for nn =1:size(datmat1,1)
        tmp = conv(datmat1(nn,:), gaussFilter);
        tmp = tmp(:,numsamps/2:end-numsamps/2);
        datmat1(nn,:) = tmp;
        
    end
    
    % ---------------- bin activity
    TrimDown = 1;
    [X, xtimes] = lt_neural_v2_QUICK_binFR(datmat1, t, binsize, TrimDown);
    
    plot(t, abs(DATMAT(1,:)));
    plot(t, datmat1(1,:), 'k', 'LineWidth', 2);
    plot(xtimes, X(1,:), 'o-r', 'LineWidth', 2);
end