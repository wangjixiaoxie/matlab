function lt_neural_commonNoise(batchf, ChansToPlot, plotOn, save_raw_dat, ...
    PlotAllChanSep, ChanToSave, onlyPlotFinalRegressionFig, skipPCAmethod, Subsample_NumSec)

%%

RegressionType = 'OLS';


%% lt 8/17/16 - exploratory only - concats neural (all chans) + song
% does not save anything
%
% ChansToPlot=neural chans (chip chan)
%% put names of all files in batch into a cell array

filenames=lt_batchsong_NamesToCell(batchf);


if onlyPlotFinalRegressionFig ==1
    plotOn=0;
    PlotAllChanSep=0;
end



%%  concatenate channels as vectors - save metadata information along the way

fs_all=[];
ampDat_all=cell(length(ChansToPlot),1); % chan x data
songDat_all=[];
transsamps=[];

% ----- COLLECT DATA
for i=1:length(filenames)
    
    fname=filenames{i};
    
    [amplifier_data,~ ,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(fname);
    
    % -- collect song
    songDat_all=[songDat_all board_adc_data(1,:)];
    transsamps=[transsamps length(board_adc_data(1,:))];
    
    % --- collect all amp dat
    for j=1:length(ChansToPlot)
        chan=ChansToPlot(j);
        ind=[amplifier_channels.chip_channel]==chan;
        dattmp=amplifier_data(ind, :); % dat for this chan
        
        dat = lt_neural_filter(dattmp, frequency_parameters); % filter bandpass
        ampDat_all{j}=[ampDat_all{j} dat]; % concat
        
        % --- collect all fs
    end
    fs_all=[fs_all frequency_parameters.amplifier_sample_rate];
    
    
end


% --- confirm that all have same sample rate
if length(fs_all)>1
    assert(all(diff(fs_all)==0), 'problem, FS not same for all file')
end

%% PLOT

if plotOn ==1
    figcount=1;
    subplotrows=length(ChansToPlot)+1;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots=[];
    
    % - song, raw
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    tt=[1:length(songDat_all)]/fs_all(1);
    plot(tt, songDat_all);
    
    %  line for each transition
    transsamps=cumsum(transsamps);
    for i=1:length(transsamps)
        x=tt(transsamps(i));
        line([x x], ylim)
    end
    
    % - plot each chan
    for i=1:length(ChansToPlot)
        chan=ChansToPlot(i);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        title(['chan ' num2str(chan)]);
        
        % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
        plot(tt, ampDat_all{i}, 'k');
    end
    
    
    linkaxes(hsplots, 'x');
end

%% ===== PCA METHOD
if skipPCAmethod==0
    % --- 1) filter each chan and putting into PCA);
    DatMat = nan(size(ampDat_all{1},2), length(ampDat_all));
    
    for i=1:length(ampDat_all)
        %    dat = lt_neural_filter(ampDat_all{i}, frequency_parameters);
        dat = ampDat_all{i};
        DatMat(:, i) = dat;
    end
    assert(~any(any(isnan(DatMat))), 'problem');
    
    % --- PCA
    [coeff, score, latent] = pca(DatMat);
    lt_figure; hold on;
    
    % -- 1) plot first few PC
    lt_subplot(4, 2, 1); hold on;
    title('coeffs (PC1)');
    xlabel('original chan');
    
    y = coeff(:,1);
    lt_plot_bar(1:length(y), y);
    
    %% -- 2) Recreate all channels only using first N PC
    if plotOn==1
        N = [1]; % which dim?
        Artifact_PCA = coeff(:, N) * score(:, N)';
        
        figcount=1;
        subplotrows=length(ChansToPlot)+1;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots=[];
        
        % - song, raw
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        tt=[1:length(songDat_all)]/fs_all(1);
        plot(tt, songDat_all);
        
        %  line for each transition
        for i=1:length(transsamps)
            x=tt(transsamps(i));
            line([x x], ylim)
        end
        
        % - plot each chan
        for i=1:length(ChansToPlot)
            chan=ChansToPlot(i);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots=[hsplots hsplot];
            title(['chan ' num2str(chan)]);
            
            % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
            dat = ampDat_all{i};
            plot(tt, dat, 'k');
            
            % -- plot reduced dimension output
            plot(tt, Artifact_PCA(i, :), 'r-');
            
            % -- plot after subtracting 1st dim
            % plot(tt, dat - tmp(i, :), 'b-');
            
            
        end
        
        
        linkaxes(hsplots, 'x');
    end
    
    %% --- PLOT EACH CHANNEL AFTER SUBTRACTING FIRST PC
    if plotOn==1;
        figcount=1;
        subplotrows=length(ChansToPlot)+1;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots=[];
        
        % - song, raw
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        tt=[1:length(songDat_all)]/fs_all(1);
        plot(tt, songDat_all);
        
        %  line for each transition
        for i=1:length(transsamps)
            x=tt(transsamps(i));
            line([x x], ylim)
        end
        
        % - plot each chan
        for i=1:length(ChansToPlot)
            chan=ChansToPlot(i);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots=[hsplots hsplot];
            title(['chan ' num2str(chan)]);
            
            % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
            dat = ampDat_all{i};
            plot(tt, dat, 'Color', [0.8 0.8 0.8]);
            
            % -- plot after subtracting 1st dim
            plot(tt, dat - Artifact_PCA(i, :), 'k-');
            
        end
        
        
        linkaxes(hsplots, 'x');
    end
end

%% ===== METHOD 2 - use other channels to predict a given channel's activity.
% NOTE: THIS IS SINGLE ELECTRODE - BELOW RUNS FOR ALL ELECTRODES.
if (0);
    ChanToExtract = 8;
    OtherChans = ChansToPlot(ChansToPlot ~= ChanToExtract);
    
    neurDat_ChanToExtract = ampDat_all{ChansToPlot == ChanToExtract}';
    
    neurDat_OtherChans = nan(size(ampDat_all{1},2), length(OtherChans));
    for i=1:length(OtherChans)
        chan = OtherChans(i);
        neurDat_OtherChans(:, i) =  ampDat_all{ChansToPlot == chan};
    end
    neurDat_OtherChans = [ones(size(neurDat_OtherChans,1),1) neurDat_OtherChans];
    
    % --- regress
    b = regress(neurDat_ChanToExtract, neurDat_OtherChans);
    
    % --- PLOT, COMPARE ACTUAL TO PREDICTED
    lt_figure; hold on;
    hsplots = [];
    
    hsplot = lt_subplot(3,1,1); hold on;
    title('song');
    hsplots=[hsplots hsplot];
    plot(tt, songDat_all);
    
    hsplot = lt_subplot(3,1,2); hold on;
    hsplots=[hsplots hsplot];
    title('PCA method');
    plot(tt, neurDat_ChanToExtract); % actual
    plot(tt, Artifact_PCA(ChansToPlot == ChanToExtract, :), 'r'); % artifact
    
    hsplot = lt_subplot(3,1,3); hold on;
    hsplots=[hsplots hsplot];
    title('regression method');
    plot(tt, neurDat_ChanToExtract); % actual
    
    Artifact_regress = b' * neurDat_OtherChans';
    plot(tt, Artifact_regress, 'r'); % artifact
    
    
    linkaxes(hsplots, 'x');
    
    
    % ---- PLOT, SUBTRACTING ARTIFACT
    lt_figure; hold on;
    hsplots = [];
    
    hsplot = lt_subplot(3,1,1); hold on;
    title('song');
    hsplots=[hsplots hsplot];
    plot(tt, songDat_all);
    
    hsplot = lt_subplot(3,1,2); hold on;
    hsplots=[hsplots hsplot];
    title('PCA method');
    plot(tt, neurDat_ChanToExtract, 'Color', [0.8 0.8 0.8]); % actual
    plot(tt, neurDat_ChanToExtract' - Artifact_PCA(ChansToPlot == ChanToExtract, :), 'r'); % artifact
    
    hsplot = lt_subplot(3,1,3); hold on;
    hsplots=[hsplots hsplot];
    title('regression method');
    plot(tt, neurDat_ChanToExtract, 'Color', [0.8 0.8 0.8]); % actual
    
    Artifact_regress = b' * neurDat_OtherChans';
    plot(tt, neurDat_ChanToExtract' - Artifact_regress, 'r'); % artifact
    
    
    linkaxes(hsplots, 'x');
end


%%  redo, but once for each channel, and collect output artifacts
Artifact_regress_all = [];
for jjj=1:length(ChansToPlot)
    ChanToExtract = ChansToPlot(jjj);
    OtherChans = ChansToPlot(ChansToPlot ~= ChanToExtract);
    
    neurDat_ChanToExtract = ampDat_all{ChansToPlot == ChanToExtract}';
    
    neurDat_OtherChans = nan(size(ampDat_all{1},2), length(OtherChans));
    for i=1:length(OtherChans)
        chan = OtherChans(i);
        neurDat_OtherChans(:, i) =  ampDat_all{ChansToPlot == chan};
    end
    neurDat_OtherChans = [ones(size(neurDat_OtherChans,1),1) neurDat_OtherChans];
    
    % --- only want weights to take into account timebins when all
    % predictor channels are fluctuation together, so penalize times when
    % one of the channels is spiking (i.e. should ignore those). do that by
    % adding a feature to the regression,
    
    % method 1) which takes into account the
    % spread of the data across electrodes at each timepoint.
    tmp = neurDat_OtherChans(:, 2:end);
    % first z-score each column
    meanmat = mean(tmp, 1);
    meanmat = repmat(meanmat, size(tmp,1),1);
    stdmat = std(tmp, 0, 1);
    stdmat = repmat(stdmat, size(tmp,1), 1);
    
    tmp_zscore = (tmp - meanmat)./stdmat;
    
    % method 2) use median value across electrodes (idea is that if is
    % artifact, then median value should be large (neg or pos), but if just
    % one predictor is spiking, then median value will not care
    median_zscored = median(tmp_zscore, 2);
    polynomial_val = prod(tmp_zscore, 2);
    
    
    % ----- ASIDE, PLOT MEDIAN TIMECOURSE
    if (0)
        figcount=1;
        subplotrows=length(ChansToPlot)+2;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots=[];
        
        % - song, raw
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        tt=[1:length(songDat_all)]/fs_all(1);
        plot(tt, songDat_all);
        
        % - plot each chan
        for i=1:length(OtherChans)
            chan=OtherChans(i);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots=[hsplots hsplot];
            title(['chan ' num2str(chan)]);
            
            % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
            dat = ampDat_all{i};
            plot(tt, dat, 'Color', [0.8 0.8 0.8]);
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        plot(tt, median_zscored, '-k'); % --- plot median zscore val
        title('median zscore');
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        plot(tt, polynomial_val, '-k'); % --- plot median zscore val
        title('polynomial (all mult)');
        
        linkaxes(hsplots, 'x');
    end
    
    %  STOPPED HERE ----- MEDIAN WORKS!! -  make sure 1) have at least 2
    %  chans (try min?), and 2) put into regression
    % at each timebin, get std across electrodes
    % TRY 1) MINIMUM, AND 2) POLYNOMIAL OF ALL CHANNELS
    %     neurDat_OtherChans = [neurDat_OtherChans median_zscored
    %     polynomial_val]; % output looks almost identical to just using
    %     median.
    neurDat_OtherChans = [neurDat_OtherChans median_zscored]; % settled on this based on looking at a few song files. works OK.
    %     neurDat_OtherChans = [neurDat_OtherChans median_zscored ...
    %         repmat(median_zscored, 1, size(neurDat_OtherChans, 2)).*neurDat_OtherChans]; % also models interaction between median and neural activity - checked and output looks almost identical to just using neural + median
    
    
    
    
    % method 2) local correlation
    
    
    % =================== IF DATAFILE TOO LONG, THEN SUBSAMPLE BEFORE
    % PERFORM REGRESSION
    % cut up into segments, get N segments, which add up to numsamples
    fs = frequency_parameters.amplifier_sample_rate;
    numsamples = Subsample_NumSec * fs;
    
        neurDat_ChanToExtract_sub = []; % these will be used to build model
        neurDat_OtherChans_sub = [];
    if length(neurDat_ChanToExtract)>numsamples
        disp('SUBSAMPLING!!');
        
        tmp = max([4, Subsample_NumSec/10]); % each chunk at minimum 4s, but try to make it 1/10 of total desired dur
        numsamples_onesegment = tmp*fs; % each segment 10 s long
        numsegmentsneeded = Subsample_NumSec/tmp;
        
        inds = 1:numsamples_onesegment:length(neurDat_ChanToExtract);
        indsToUse = 1:round(length(inds)/numsegmentsneeded):length(inds);
        
        inds = inds(indsToUse);
        
        for j=1:length(inds)
            
            indstmp = inds(j):inds(j)+numsamples_onesegment; % sample inds to take
            
            if max(indstmp)>length(neurDat_ChanToExtract)
                neurDat_ChanToExtract_sub = [neurDat_ChanToExtract_sub; ...
                    neurDat_ChanToExtract(indstmp(1):end)];
                neurDat_OtherChans_sub = [neurDat_OtherChans_sub; ...
                    neurDat_OtherChans(indstmp(1):end, :)];
            else
                neurDat_ChanToExtract_sub = [neurDat_ChanToExtract_sub; ...
                    neurDat_ChanToExtract(indstmp)];
                neurDat_OtherChans_sub = [neurDat_OtherChans_sub; ...
                    neurDat_OtherChans(indstmp, :)];
            end
            
            disp([indstmp(1) indstmp(end)]./fs);
        end
        
        
    else
        % then use all data.
        neurDat_ChanToExtract_sub = neurDat_ChanToExtract;
        neurDat_OtherChans_sub = neurDat_OtherChans;
    end
    
    % --- regress
    if strcmp(RegressionType, 'OLS');
        b = regress(neurDat_ChanToExtract_sub, neurDat_OtherChans_sub);
    elseif strcmp(RegressionType, 'ridge')
        % - ridge regression
        % D = x2fx(neurDat_OtherChans,'interaction');
        % D(:,1) = []; % No constant term
        k = [0 1e-4 5e-3];
        b = ridge(neurDat_ChanToExtract_sub, neurDat_OtherChans_sub,k);
        % - plot ridge
        figure
        plot(k,b,'LineWidth',2)
        ylim([-100 100])
        grid on
        xlabel('Ridge Parameter')
        ylabel('Standardized Coefficient')
        title('{\bf Ridge Trace}')
        legend('x1','x2','x3','x1x2','x1x3','x2x3')
        b = b(:,2);
    elseif strcmp(RegressionType, 'lasso')
        
    end
    
    
    
    
    
    % --- extract artifact
    Artifact_regress = b' * neurDat_OtherChans';
    Artifact_regress_all = [Artifact_regress_all; Artifact_regress];
end

%% === PLOT ALL CHANNELS, regression method, overlaying artifact.

if plotOn==1
    % ================ METHOD 2 (regression)
    figcount=1;
    subplotrows=length(ChansToPlot)+1;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots=[];
    
    % - song, raw
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    tt=[1:length(songDat_all)]/fs_all(1);
    plot(tt, songDat_all);
    
    %  line for each transition
    for i=1:length(transsamps)
        x=tt(transsamps(i));
        line([x x], ylim)
    end
    
    % - plot each chan
    for i=1:length(ChansToPlot)
        chan=ChansToPlot(i);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        title(['chan ' num2str(chan)]);
        
        % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
        dat = ampDat_all{i};
        plot(tt, dat, 'Color' , 'k');
        
        % -- plot after subtracting 1st dim
        plot(tt, Artifact_regress_all(i, :), 'r-');
        
    end
    
    
    linkaxes(hsplots, 'x');
    lt_subtitle('[data + artifact] regression method');
    
end

%% ====== PLOT EACH CHANNEL, COMPARING BOTH METHODS.
% NOTE: CUTS DATA SHORT IF TOO LONG (GREATER THAN 12000000 SAMPS) SO SAVED
% FILE IS NOT TOO LARGE

if plotOn ==1 || onlyPlotFinalRegressionFig ==1
    % ================ METHOD 1 (PCA)
    if skipPCAmethod==0
        figcount=1;
        subplotrows=length(ChansToPlot)+1;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots=[];
        
        % - song, raw
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        tt=[1:length(songDat_all)]/fs_all(1);
        plot(tt, songDat_all);
        
        %  line for each transition
        for i=1:length(transsamps)
            x=tt(transsamps(i));
            line([x x], ylim)
        end
        
        % - plot each chan
        for i=1:length(ChansToPlot)
            chan=ChansToPlot(i);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots=[hsplots hsplot];
            title(['chan ' num2str(chan)]);
            
            % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
            dat = ampDat_all{i};
        if length(dat) > 12000000
            dat = dat(1:12000000);
        plot(tt, dat, 'Color', [0.3 0.8 0.8]);
        end
            
            % -- plot after subtracting 1st dim
            plot(tt, dat - Artifact_PCA(i, 1:length(dat)), 'k-');
            
        end
        
        lt_subtitle('PCA method [final + original]');
        linkaxes(hsplots, 'x');
    end
    
    % ================ METHOD 2 (regression)
    figcount=1;
    subplotrows=length(ChansToPlot)+1;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots=[];
    
    % - song, raw
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots=[hsplots hsplot];
    tt=[1:length(songDat_all)]/fs_all(1);
    plot(tt, songDat_all);
    
    %  line for each transition
    for i=1:length(transsamps)
        x=tt(transsamps(i));
        line([x x], ylim)
    end
    
    % - plot each chan
    for i=1:length(ChansToPlot)
        chan=ChansToPlot(i);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        title(['chan ' num2str(chan)]);
        
        % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
        dat = ampDat_all{i};
        if length(dat) > 12000000
            dat = dat(1:12000000);
            tt = tt(1:12000000);
        end
        plot(tt, dat, 'Color', [0.3 0.8 0.8]);
        
        % -- plot after subtracting 1st dim
        plot(tt, dat - Artifact_regress_all(i, 1:length(dat)), 'k-');
    end
    
    
    linkaxes(hsplots, 'x');
    lt_subtitle('regression method [final + original]');
    
end

%% ====== ONE PLOT FOR EACH CHANNEL, ONE SUBPLOT FOR EACH METHOD.
% NOTE: to do this stop here and by hand plot each one separately by
% changing value of i

if PlotAllChanSep ==1
    
    for i=1:length(ChansToPlot)
        
        chan=ChansToPlot(i);
        figcount=1;
        subplotrows=3;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots=[];
        
        % 1) - song, raw
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        tt=[1:length(songDat_all)]/fs_all(1);
        plot(tt, songDat_all);
        
        %  line for each transition
        for ii=1:length(transsamps)
            x=tt(transsamps(ii));
            line([x x], ylim)
        end
        
        
        % 2) ================ METHOD 1 (PCA)
        if (1)
            % - plot each chan
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots=[hsplots hsplot];
            title(['PCA (orig + new), chan ' num2str(chan)]);
            
            % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
            dat = ampDat_all{i};
            %         plot(tt, dat, 'Color', [0.8 0.8 0.8]);
            
            % -- plot after subtracting 1st dim
            plot(tt, dat, 'Color', [0.3 0.8 0.8]);
            plot(tt, dat - Artifact_PCA(i, :), 'k-');
        end
        
        % 3) =============== METHOD 2, Regression
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots=[hsplots hsplot];
        title(['regression (orig + new), chan ' num2str(chan)]);
        
        % dat=lt_neural_filter(ampDat_all{i}, frequency_parameters);
        dat = ampDat_all{i};
        %         plot(tt, dat, 'Color', [0.8 0.8 0.8]);
        
        % -- plot after subtracting 1st dim
        plot(tt, dat, 'Color', [0.3 0.8 0.8]);
        plot(tt, dat - Artifact_regress_all(i, :), 'k-');
        
        linkaxes(hsplots, 'x');
        
    end
end


%% ===== SAVE A SPECIFIC CHANNEL, TO USE FOR CLUSTERING
% USES REGRESSION OUTPUT.

if save_raw_dat == 1
    
    % === COLLECT METADAT
    metaDat=struct;
    % ----- COLLECT DATA
    for i=1:length(filenames)
        
        fname=filenames{i};
        
        [frequency_parameters, board_adc_data] = pj_readIntanNoGui_AudioOnly(fname);
        
        % --- collect metadata
        metaDat(i).filename=fname; % filename
        metaDat(i).numSamps=length(board_adc_data(1,:)); % length of file (samps)
        metaDat(i).fs=frequency_parameters.amplifier_sample_rate; % fs
        metaDat(i).songDat=board_adc_data(1,:);
    end
    
    % --- save a data file with name corresponding
    dirsavename=['Chan' num2str(ChanToSave) 'amp-' batchf];
    
    currdir=pwd;
    
    % --- extract data
    ind = ChansToPlot == ChanToSave;
    data = ampDat_all{ind}; % wave_clus needs this name change
    data = data - Artifact_regress_all(ind, :);
    
    % --- make a new dir for this data
    if ~exist(dirsavename, 'dir')
        % make dir
        mkdir(dirsavename)
        cd(dirsavename)
        
        
        % --- save
        save(['data.mat'], 'data', '-v7.3');  % save vector
        save(['MetaDat.mat'], 'metaDat', '-v7.3'); % save metadata
        eval(['!cp ../' batchf ' .']) % copy batchfile over
    else
        if strcmp(input('data already exists - overwrite folder (y) or just data.mat (n)?' , 's'), 'y')
            eval(['!rm -r ' dirsavename ]);
            mkdir(dirsavename)
            cd(dirsavename)
            
            disp(['... deleted folder ' dirsavename]);
            
            % --- save
            save(['data.mat'], 'data', '-v7.3');  % save vector
            save(['MetaDat.mat'], 'metaDat', '-v7.3'); % save metadata
            eval(['!cp ../' batchf ' .']) % copy batchfile over
        else
            % then save data (overwrite)
            cd(dirsavename)
            save(['data.mat'], 'data', '-v7.3');  % save vector
            %     save(['MetaDat.mat'], 'metaDat'); % save metadata
        end
    end
    
    tstamp = lt_get_timestamp(0);
    fid = fopen(['DENOISED_' tstamp], 'w');
    fclose(fid);
    
    % --- SAVE FIGS
    try cd('FIGS')
    catch err
        mkdir('FIGS')
    end
    cd('FIGS');
    lt_save_all_figs;
    
    
    cd(currdir)
    
end





