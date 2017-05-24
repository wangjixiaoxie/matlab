%% LT 8/30/16 - given segments, 1) calcualtes mean firing rate and other stats in window pre and post and 2) outputs smoothed firing rate
function [FiringRateOut BurstFracOut]=lt_neural_CalMeanFiring(SegmentsExtract, Params, clustwanted, Window_relOnset, motifwind, Window_relOffset, ISIthreshburst)
% outputs stats for:
% 2) entire motif [DEFAULT]

% REQUIRED:
% motifwind=[-0.04 -0.04]; % e.g. 40ms premotif onset to 40ms pre motif offest

% OPTIONAL: [make sure these are in temporal order]
% Window_relOnset{1}=[-20 -18]; % window rel 1st syl onset
% Window_relOnset{2}=[-3 -2.5]; % window rel 1st syl onset
% Window_relOnset{3}=[-0.6 -0.1]; % and so on, can have as many windows as desired
%
% Window_relOffset{1}=[0.1 0.6]; % window rel last syl offset... can have as many as desired
% Window_relOffset{2}=[2.5 3];
% Window_relOffset{3}=[18 20];
%
% clustwanted=1; % wave_clus cluster

% --- OUTPUTS

% running average
% firing rates

% ==== PARAMS
useRescaled=0; % use linearly warped data or no.
% window_sm=0.02;
% windshift=0.002;

%%


if useRescaled==1
    spktimefield='spk_Times_scaled';
else
    spktimefield='spk_Times';
end

%% extract spike times for each trial
numtrials=length(SegmentsExtract);
Yspks_entiresegment={}; % one ind for each seg.

WindowedDat=struct; % collect spk times that are within the desired pre and post windows

for kk=1:numtrials
    inds=find([SegmentsExtract(kk).spk_Clust]==clustwanted);
    spktimes=SegmentsExtract(kk).(spktimefield)(inds);
    Yspks_entiresegment{kk}=spktimes;
    
    % -- FIND ONSET WINDOWED SPIKES
    for j=1:length(Window_relOnset)
        wind=Window_relOnset{j};
        wind=wind + Params.REGEXP.predur; % find times relative to inds, by adding predur
        spktimes_wind= spktimes(spktimes>wind(1) & spktimes<wind(2));
        assert(Params.REGEXP.predur>=-Window_relOnset{j}(1), 'PROBLEM - not enough post dur data to get wind'); % maek sure have enough data to do this
        
        % output
        WindowedDat.relMotifOnset(j).windTimes_thistrial(kk, :)=wind;
        WindowedDat.relMotifOnset(j).spktimes_inWind{kk}=spktimes_wind;
        WindowedDat.relMotifOnset(j).numspikes(kk)=length(spktimes_wind);
        WindowedDat.relMotifOnset(j).Window_relOnset=Window_relOnset{j};
    end
    
    
    % -- FIND OFFSET WINDOWED SPIKES
    % for this trial (kk), what is the time of offset of last syl in motif?
    postdur=Params.REGEXP.postdur;
    timeOfMotifOffset=SegmentsExtract(kk).global_offtime_motifInclFlank - ...
        SegmentsExtract(kk).global_ontime_motifInclFlank - postdur; % time of offset of last syl in motif
    
    % make sure have enough data to do this
    for j=1:length(Window_relOffset)
        wind=Window_relOffset{j};
        wind=wind+timeOfMotifOffset;
        assert(postdur>=Window_relOffset{j}(2), 'PROBLEM - not enough post dur data to get wind'); % maek sure have enough data to do this
        
        spktimes_wind= spktimes(spktimes>wind(1) & spktimes<wind(2));
        
        % optput
        WindowedDat.relMotifOffset(j).windTimes_thistrial(kk, :)=wind;
        WindowedDat.relMotifOffset(j).spktimes_inWind{kk}=spktimes_wind;
        WindowedDat.relMotifOffset(j).numspikes(kk)=length(spktimes_wind);
        WindowedDat.relMotifOffset(j).Window_relOffset=Window_relOffset{j};
    end
    
    
    % -- FIND SPIKES DURING MOTIF
    wind(1)=Params.REGEXP.predur+motifwind(1);
    wind(2)=timeOfMotifOffset+motifwind(2);
    
    spktimes_wind= spktimes(spktimes>wind(1) & spktimes<wind(2));
    % optput
    WindowedDat.entireMotif.windTimes_thistrial(kk, :)=wind;
    WindowedDat.entireMotif.spktimes_inWind{kk}=spktimes_wind;
    WindowedDat.entireMotif.numspikes(kk)=length(spktimes_wind);
    WindowedDat.entireMotif.Window_relOnsetOffset=motifwind;
end

%% ======== DIAGNOSTIC
if (0)
    % THIS IS USEFUL, PLOTS EACH TRIAL (SET i) SONG, SEGMENTS EXTRACTED,
    % AND SPIKES
    
    % for each trial, plot raster, plus data for segments
    lt_figure; hold on;
    i=15;
    
    Yspks_test={};
    
    % ---- PLOT FOR THIS TRIAL
    spktimes=SegmentsExtract(i).spk_Times(SegmentsExtract(i).spk_Clust == clustwanted);
    
    % -- plot spikes
    Yspks_test{1}=spktimes;
    [xbin, ymean, ysem, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks_test, 0.15, 0.05, 1, 'k');
    plot(xbin, ymean_hz, '-r');
    
    % -- plot song segment
    songdat=SegmentsExtract(i).songdat;
    tt=(1:length(songdat))/SegmentsExtract(i).fs;
    plot(tt, songdat, 'y')
    
    % -- plot motif
    DATSTRUCT=WindowedDat.entireMotif;
    stimes=DATSTRUCT.spktimes_inWind{i};
    plot(stimes, 3, 'ob');
    spnum=DATSTRUCT.numspikes(i);
    trialon=DATSTRUCT.windTimes_thistrial(i, 1);
    trialoff=DATSTRUCT.windTimes_thistrial(i, 2);
    dur=trialoff - trialon;
    srate=spnum/dur;
    line([trialoff trialon], [srate srate]);
    line([trialoff trialoff], ylim, 'Color', 'm');
    line([trialon trialon], ylim, 'Color', 'm');
    
    % - plot onset
    DATSTRUCT=WindowedDat.relMotifOnset(1);
    stimes=DATSTRUCT.spktimes_inWind{i};
    plot(stimes, 3, 'ob');
    spnum=DATSTRUCT.numspikes(i);
    trialon=DATSTRUCT.windTimes_thistrial(i, 1);
    trialoff=DATSTRUCT.windTimes_thistrial(i, 2);
    dur=trialoff - trialon;
    srate=spnum/dur;
    line([trialoff trialon], [srate srate]);
    line([trialoff trialoff], ylim, 'Color', 'm');
    line([trialon trialon], ylim, 'Color', 'm');
    
    % - plot offset
    DATSTRUCT=WindowedDat.relMotifOffset(1);
    stimes=DATSTRUCT.spktimes_inWind{i};
    plot(stimes, 3, 'ob');
    spnum=DATSTRUCT.numspikes(i);
    trialon=DATSTRUCT.windTimes_thistrial(i, 1);
    trialoff=DATSTRUCT.windTimes_thistrial(i, 2);
    dur=trialoff - trialon;
    srate=spnum/dur;
    line([trialoff trialon], [srate srate]);
    line([trialoff trialoff], ylim, 'Color', 'm');
    line([trialon trialon], ylim, 'Color', 'm');
end

%% ==== combine all trials, plot pre, motif, and post (firing rate)

Xall_times=[]; % trials x windows
Yall_rates=[]; % trials x windows
Xall_meantimes=[]; % each ind is one segment - put them in temporal order
% Yall_mean=[];
% Yall_sem=[];


% === for each WINDOWED SEGMENT
for j=1:length(WindowedDat.relMotifOnset)
    DATSTRUCT=WindowedDat.relMotifOnset(j);
    
    % -- collect
    alldurs=DATSTRUCT.windTimes_thistrial(:,2) ...
        - DATSTRUCT.windTimes_thistrial(:,1);
    allsprates=DATSTRUCT.numspikes'./alldurs;
    allmeantimes=mean(DATSTRUCT.windTimes_thistrial, 2);
    
    Xall_meantimes = [Xall_meantimes mean(allmeantimes)];
    Yall_rates  = [Yall_rates allsprates];
    Xall_times = [Xall_times allmeantimes];
end

% ==== DUR MOTIF
DATSTRUCT=WindowedDat.entireMotif;

% -- collect
alldurs=DATSTRUCT.windTimes_thistrial(:,2) ...
    - DATSTRUCT.windTimes_thistrial(:,1);
allsprates=DATSTRUCT.numspikes'./alldurs;
allmeantimes=mean(DATSTRUCT.windTimes_thistrial, 2);

Xall_meantimes = [Xall_meantimes mean(allmeantimes)];
Yall_rates  = [Yall_rates allsprates];
Xall_times = [Xall_times allmeantimes];


% ===== POST MOTIF
for j=1:length(WindowedDat.relMotifOffset)
    
    DATSTRUCT=WindowedDat.relMotifOffset(j);
    
    % -- collect
    alldurs=DATSTRUCT.windTimes_thistrial(:,2) ...
        - DATSTRUCT.windTimes_thistrial(:,1);
    allsprates=DATSTRUCT.numspikes'./alldurs;
    allmeantimes=mean(DATSTRUCT.windTimes_thistrial, 2);
    
    Xall_meantimes = [Xall_meantimes mean(allmeantimes)];
    Yall_rates  = [Yall_rates allsprates];
    Xall_times = [Xall_times allmeantimes];
end

% +++++++++++++++++ PLOT
x=size(Xall_times, 2);
lt_figure; hold on;

% ==== all trials (using actual timepoint)
lt_subplot(3,1,1); hold on;
title('firing rate, vs. mean time of window');
plot(Xall_times', Yall_rates', '-', 'Color', [0.5 0.5 0.5]);
% - mean
xmean=mean(Xall_times);
ymean=mean(Yall_rates);
ysem=lt_sem(Yall_rates);
lt_plot(xmean, ymean, {'Errors', ysem, 'LineStyle', '-'})


% ==== all trials (not real times)
lt_subplot(3,1,2); hold on;
title('firing rate, each window');
% each trial
plot(1:x, Yall_rates', '-', 'Color', [0.5 0.5 0.5]);
% - mean
xmean=mean(Xall_times);
ymean=mean(Yall_rates);
ysem=lt_sem(Yall_rates);
lt_plot(1:x, ymean, {'Errors', ysem, 'LineStyle', '-'})


% === summary
lt_subplot(3,1,3); hold on;
title('firing rate');
boxplot(Yall_rates)

% ======== for output
FiringRateOut.xmean_sec=xmean;
FiringRateOut.Yall_rates=Yall_rates;


%% ================= ISI DISTRIBUTIONS

figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

tmp=1.2;
xtmp=tmp.^(0:100);
xcenters=xtmp(xtmp<1000);

Yall_FracSpInbursts=[];


% +++++++++++++++++++++++++++++++++++++++++++++++ HISTOGRAMS
% === for each WINDOWED SEGMENT
windowname='relMotifOnset';
for j=1:length(WindowedDat.(windowname))
    
    DATSTRUCT=WindowedDat.(windowname)(j);
    
    % -- collect
    tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
    allISI=cell2mat(tmp);
    % -- fraction of neurons in bursts
    fracbursts=sum(allISI<ISIthreshburst)/length(allISI);
    Yall_FracSpInbursts=[Yall_FracSpInbursts fracbursts];
    % -- plot
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['segment: ' windowname]);
    allISI=allISI*1000; % convert to ms
    lt_plot_histogram(allISI, xcenters, 1, 1, 1, 1, 'k');
    logx
    xlim([0 1000]);
    set(gca, 'XTick', [1 5 10 100 1000]);
end

% ==== DUR MOTIF
windowname='entireMotif';
DATSTRUCT=WindowedDat.(windowname);

% -- collect
tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
allISI=cell2mat(tmp);
% -- fraction of neurons in bursts
fracbursts=sum(allISI<ISIthreshburst)/length(allISI);
Yall_FracSpInbursts=[Yall_FracSpInbursts fracbursts];
% -- plot
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['segment: ' windowname]);
allISI=allISI*1000; % convert to ms
lt_plot_histogram(allISI, xcenters, 1, 1, 1, 1, 'k');
logx
xlim([0 1000]);
set(gca, 'XTick', [1 5 10 100 1000]);


% ===== POST MOTIF
windowname='relMotifOffset';
for j=1:length(WindowedDat.(windowname))
    
    DATSTRUCT=WindowedDat.(windowname)(j);
    
    % -- collect
    tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
    allISI=cell2mat(tmp);
    % -- fraction of neurons in bursts
    fracbursts=sum(allISI<ISIthreshburst)/length(allISI);
    Yall_FracSpInbursts=[Yall_FracSpInbursts fracbursts];
    % -- plot
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['segment: ' windowname]);
    allISI=allISI*1000; % convert to ms
    lt_plot_histogram(allISI, xcenters, 1, 1, 1, 1, 'k');
    logx
    xlim([0 1000]);
    set(gca, 'XTick', [1 5 10 100 1000]);
end



% +++++++++++++++++++++++++++++++++ COMBINE ALL INTO CDF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all into cdf (num put in random y position)');
numwinds=10;
plotcols=lt_make_plot_colors(numwinds, 0, 0);
ccount=1;

% === for each WINDOWED SEGMENT
windowname='relMotifOnset';
for j=1:length(WindowedDat.(windowname))
    
    DATSTRUCT=WindowedDat.(windowname)(j);
    
    % -- collect
    tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
    allISI=cell2mat(tmp);
    % -- plot
    allISI=allISI*1000; % convert to ms
    lt_plot_cdf(allISI, plotcols{ccount}, 0);
    lt_plot_text(50+600*rand, 0.2 , num2str(ccount), plotcols{ccount});
    logx; xlim([0 1000]);
    set(gca, 'XTick', [1 5 10 100 1000]);
    
    ccount=ccount+1;
end

% ==== DUR MOTIF
windowname='entireMotif';
DATSTRUCT=WindowedDat.(windowname);

% -- collect
tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
allISI=cell2mat(tmp);
% -- plot
allISI=allISI*1000; % convert to ms
lt_plot_cdf(allISI, plotcols{ccount}, 0);
lt_plot_text(50+600*rand, 0.2 , num2str(ccount), plotcols{ccount});
logx; xlim([0 1000]);
set(gca, 'XTick', [1 5 10 100 1000]);

ccount=ccount+1;


% ===== POST MOTIF
windowname='relMotifOffset';
for j=1:length(WindowedDat.(windowname))
    
    DATSTRUCT=WindowedDat.(windowname)(j);
    
    % -- collect
    tmp=cellfun(@diff, DATSTRUCT.spktimes_inWind, 'UniformOutput', 0); % get all isi
    allISI=cell2mat(tmp);
    % -- plot
    allISI=allISI*1000; % convert to ms
    lt_plot_cdf(allISI, plotcols{ccount}, 0);
    lt_plot_text(50+600*rand, 0.2 , num2str(ccount), plotcols{ccount});
    logx; xlim([0 1000]);
    set(gca, 'XTick', [1 5 10 100 1000]);
    
    ccount=ccount+1;
end


% ======================== PLOT FRAC OF NEUROSNIN BURSTS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('frac spikes in bursts');
xlabel('window');
plot(Yall_FracSpInbursts, '-ok');
ylim([0 1]);
BurstFracOut.Yall_FracSpInbursts=Yall_FracSpInbursts;


%% -
% [xbin, ymean, ysem, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks_entiresegment, window_sm, windshift, 0, '');



