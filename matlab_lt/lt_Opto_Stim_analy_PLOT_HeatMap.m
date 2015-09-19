%% LT 3/9/15 - compare two conditions (e.g. stim vs. stim catch) that both have stim timing, looking at effects triggered by stim.

function [Params, DATSTRUCT]=lt_Opto_Stim_analy_PLOT_HeatMap(Params,StatsStruct);

fieldname1='StimCatch'; % used as "control"
fieldname2='StimNotCatch'; % test group


%% AUTO PARAMS
% conversion between PC samples and time (ms)

% Extract old params
NumFields=length(Params.FieldsToCheck);
PC_freqrange=Params.PC_freqrange;
pc_harms=Params.pc_harms;
Fs=Params.Fs;
PreDur=1000*Params.PreDur; % time of alignment (i.e. duration of data taken before target syl), ms.
StimDur=Params.StimDur;


%% 1) get median or mean PC for the non-stim trials


fieldname=fieldname1;

PC=StatsStruct.(fieldname).PC;

DATSTRUCT.(fieldname).PC_mean=mean(PC,2);
DATSTRUCT.(fieldname).PC_median=median(PC,2);

% for each time bin, get std (i.e. to get z-score later);
DATSTRUCT.(fieldname).PC_std=std(PC,0,2);




%% 2) for all stim trials, get deviation from the non-stim, for each time bin

fieldname=fieldname2;

% FIRST, get all stim trials into a matrix, where each row is a trial
DATSTRUCT.(fieldname).PC=StatsStruct.(fieldname).PC';

% for each trial, subtract median for non-stim trials
NumStimTrials=size(DATSTRUCT.(fieldname).PC,1);

DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC=DATSTRUCT.(fieldname).PC - repmat(DATSTRUCT.StimCatch.PC_mean',NumStimTrials,1);  % mean
DATSTRUCT.(fieldname).PC_SubtractNonstimMednPC=DATSTRUCT.(fieldname).PC- repmat(DATSTRUCT.StimCatch.PC_median',NumStimTrials,1);% subtract median

% then divide by std of each time bin to get z-score
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean=DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC./repmat(DATSTRUCT.StimCatch.PC_std',NumStimTrials,1);

% Throw out data outside of epochs of interest (as high noise overshadows
% signal).
NumTimeWindows=length(Params.TimeWindowList);
% convert time windows to samples
for i=1:NumTimeWindows;
    ind1_ms=Params.TimeWindowList{i}(1); % in ms
    ind2_ms=Params.TimeWindowList{i}(2);
    
    [~,ind1_samp]=min(abs(1000*Params.tf_bins.tPC-ind1_ms)); % convert to sample num
    [~,ind2_samp]=min(abs(1000*Params.tf_bins.tPC-ind2_ms));
   
    Params.HEATMAP.TimeWindowList_tPC{i}=[ind1_samp ind2_samp];
end

% only keep those epochs
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed=nan(size(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean)); % get empty matrix of nans
DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed=nan(size(DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC)); % get empty matrix of nans
for i=1:NumTimeWindows;
    t1=Params.HEATMAP.TimeWindowList_tPC{i}(1);
    t2=Params.HEATMAP.TimeWindowList_tPC{i}(2);
    DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(:,t1:t2)=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean(:,t1:t2);
    DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed(:,t1:t2)=DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC(:,t1:t2);

end


%% 2.5) Do the same for stim catch trials
fieldname=fieldname1;

% FIRST, get all stim trials into a matrix, where each row is a trial
DATSTRUCT.(fieldname).PC=StatsStruct.(fieldname).PC';

% for each trial, subtract median for non-stim trials
NumTrials=size(DATSTRUCT.(fieldname).PC,1);

DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC=DATSTRUCT.(fieldname).PC - repmat(DATSTRUCT.StimCatch.PC_mean',NumTrials,1);  % mean
DATSTRUCT.(fieldname).PC_SubtractNonstimMednPC=DATSTRUCT.(fieldname).PC- repmat(DATSTRUCT.StimCatch.PC_median',NumTrials,1);% subtract median

% then divide by std of each time bin to get z-score
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean=DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC./repmat(DATSTRUCT.StimCatch.PC_std',NumTrials,1);

% Throw out data outside of epochs of interest (as high noise overshadows
% signal).
NumTimeWindows=length(Params.TimeWindowList);

% only keep those epochs
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed=nan(size(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean)); % get empty matrix of nans
DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed=nan(size(DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC)); % get empty matrix of nans

for i=1:NumTimeWindows;
    t1=Params.HEATMAP.TimeWindowList_tPC{i}(1);
    t2=Params.HEATMAP.TimeWindowList_tPC{i}(2);
    DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(:,t1:t2)=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean(:,t1:t2);
    DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed(:,t1:t2)=DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC(:,t1:t2);
end



%% 3) sort all stim trials by onset timing of stim (row 1 1st, last row last)
fieldname=fieldname2;

% sort by stim onset time, and get matrix of stim onset times.
[TimeSinceLastTrig_sorted, inds]=sort(StatsStruct.(fieldname).TimeSinceLastTrig);
TrigTimes_fromEpochStart=PreDur-TimeSinceLastTrig_sorted;
 
TrigTimes_fromEpochStart=TrigTimes_fromEpochStart+Params.HEATMAP.stim_delay;  % add delay between trigger and stim onset

% OUTPUT:
Params.HEATMAP.TrigTimes_fromEpochStart_delayed=TrigTimes_fromEpochStart;



% for i=1:length(TrigTimes_fromEpochStart);
%     tt=TrigTimes_fromEpochStart(i);
%     [~,ind]=min(abs(1000*Params.tf_bins.tPC-tt)); % convert trig times to sample times
% TrigTimes_fromEpochStart_samples(i)=ind;
% end

% use stim onset sorting indices to sort the PC trials.
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(inds,:);
DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed=DATSTRUCT.(fieldname).PC_SubtractNonstimMeanPC_windowed(inds,:);


%% PLOT Z-score
% concatenate [stim;nonstim]
DATSTRUCT.All.PC_ZscoreRelNonstimMean_windowed=[DATSTRUCT.StimNotCatch.PC_ZscoreRelNonstimMean_windowed; ...
    DATSTRUCT.StimCatch.PC_ZscoreRelNonstimMean_windowed];


% PLOT
% clims=[-3 3];
% figure; hold on;
% imagesc(Params.tf_bins.tPC,1:NumStimTrials,DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed,clims)
% colormap('spring');
% title('Stim trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
% xlabel('time');
% ylabel('trial #')
% zlabel('z-score')
% colorbar
% 
% for i=1:NumStimTrials;
%     plot(TrigTimes_fromEpochStart(i)/1000,i,'kx');
% end
    
N=size(DATSTRUCT.All.PC_ZscoreRelNonstimMean_windowed,1);

clims=[-3 3];
figure; hold on;
imagesc(Params.tf_bins.tPC,1:N,DATSTRUCT.All.PC_ZscoreRelNonstimMean_windowed,clims)
colormap('spring');
title('Each trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
zlabel('z-score')
colorbar

for i=1:length(TrigTimes_fromEpochStart);
    plot(TrigTimes_fromEpochStart(i)/1000,i,'kx');
end



% % ONLY PLOT INDIVIDUAL TIME WINDOWS
% t1=Params.HEATMAP.TimeWindowList_tPC{2}(1);
% t2=Params.HEATMAP.TimeWindowList_tPC{2}(2);
% 
% figure; hold on;
% imagesc(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(:,t1:t2),clims)
% colormap('spring')
% 

%% PLOT FF DIFF

DATSTRUCT.All.PC_SubtractNonstimMeanPC_windowed=[DATSTRUCT.StimNotCatch.PC_SubtractNonstimMeanPC_windowed; ...
    DATSTRUCT.StimCatch.PC_SubtractNonstimMeanPC_windowed];

C(2)=max(max(DATSTRUCT.All.PC_SubtractNonstimMeanPC_windowed));
C(1)=min(min(DATSTRUCT.All.PC_SubtractNonstimMeanPC_windowed));

clims=C;
figure; hold on;
imagesc(Params.tf_bins.tPC,1:N,DATSTRUCT.All.PC_SubtractNonstimMeanPC_windowed,clims)
colormap('spring');
title('Each trial pitch, FF diff relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
colorbar

for i=1:length(TrigTimes_fromEpochStart);
    plot(TrigTimes_fromEpochStart(i)/1000,i,'kx');
end
    



%% 4) plot heat maps, but square z-score first, (i.e. magnify differences, and negatives)

if (0)
DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed_squared=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed.^2;

% PLOT
figure; hold on;
imagesc(Params.tf_bins.tPC,1:NumStimTrials,log(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed_squared),clims)
colormap('spring');
title('Stim trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
colorbar

for i=1:NumStimTrials;
    plot(TrigTimes_fromEpochStart(i)/1000,i,'ks');
end


Y=log(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed_squared);

tPC=Params.tf_bins.tPC;

figure; hold on;
    plot(tPC*1000,Y,'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade

    plot(tPC*1000,mean(Y),'Linewidth',2)

    line(xlim,[0 0])


end

%% RUNNING AVERAGE PC (over trials, sorted by stim onset)

NumTimeBins=size(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed,2);

fieldnameList={'StimNotCatch','StimCatch'};
for kk=1:length(fieldnameList);
    fieldname=fieldnameList{kk};
for i=1:NumTimeBins; % for all time bins
    DATSTRUCT.(fieldname).trialsmoothed.PC_ZscoreRelNonstimMean_windowed(:,i)=smooth(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(:,i),Params.HEATMAP.runningbin); % smooth PCs
   
    
    Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed=smooth(Params.HEATMAP.TrigTimes_fromEpochStart_delayed, Params.HEATMAP.runningbin); % smooth trigger times
end

    % remove start and end as smooth binsize lowers (until is 1 at
    % endpoints)
trialstoremove=floor(Params.HEATMAP.runningbin/2); % how many from start and end
    
DATSTRUCT.(fieldname).trialsmoothed.PC_ZscoreRelNonstimMean_windowed(1:trialstoremove,:)=[];
DATSTRUCT.(fieldname).trialsmoothed.PC_ZscoreRelNonstimMean_windowed(end-trialstoremove+1:end,:)=[];
Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed(1:trialstoremove)=[];
Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed(end-trialstoremove+1:end)=[];
end
  
% combine stim and nostim
DATSTRUCT.All.trialsmoothed.PC_ZscoreRelNonstimMean_windowed=[DATSTRUCT.StimNotCatch.trialsmoothed.PC_ZscoreRelNonstimMean_windowed;...
    DATSTRUCT.StimCatch.trialsmoothed.PC_ZscoreRelNonstimMean_windowed];


% PLOT
N=size(DATSTRUCT.All.trialsmoothed.PC_ZscoreRelNonstimMean_windowed,1); % num trials

% clims=[-3 3];
figure; hold on;
imagesc(Params.tf_bins.tPC,1:N,DATSTRUCT.All.trialsmoothed.PC_ZscoreRelNonstimMean_windowed)
colormap('spring');
title(['All trial pitch, z-scored relative to non-stim trial mean and STD; binsize: ' num2str(Params.HEATMAP.runningbin)]);
xlabel('time');
ylabel('trial #')
colorbar

for i=1:length(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed);
        plot(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed(i)/1000,i,'kx');
end


    


%% do same, but normalizing each trial's z-score to itself (i.e. divide by max z-score)

% if (0)
% tmp=max(abs(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed(inds,:)),[],2); % gets max zscore for each trial (in time windows only);
% 
% DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed_WithinTrialNorm=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed./repmat(tmp,1,size(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed,2));
% 
% 
% 
% 
% % PLOT
% figure; hold on;
% imagesc(DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed_WithinTrialNorm)
% end

%% HYPOTHESIS TESTING

% Permutation test




%% PLOTTING TEST
fieldname=fieldname2;

Y=DATSTRUCT.(fieldname).PC_ZscoreRelNonstimMean_windowed;

tPC=Params.tf_bins.tPC;

figure; hold on;
    plot(tPC*1000,Y,'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade

    plot(tPC*1000,mean(Y),'Linewidth',2)

    line(xlim,[0 0])
title('Stim trial pitch contours, each time bin z-scored rel. to catch mean PC. All trials + mean');
ylabel('z-score');
xlabel('time (ms)');

    
    
%% SAME PLOT, but for amplitude profile - interesting point to look tempo changes

fieldname=fieldname2;
% First plot all trials, simply raw values

Params.HEATMAP.TrigTimes_fromEpochStart_delayed; % trig times are same


%% SORT STIM TRIALS BY STIM TIME
[~, inds]=sort(StatsStruct.(fieldname).TimeSinceLastTrig);


% DATSTRUCT.(fieldname2).sm_log= StatsStruct.(fieldname2).sm_log';
sm=StatsStruct.(fieldname2).sm;
sm_log=log(sm);
DATSTRUCT.(fieldname2).sm_log= sm_log';

DATSTRUCT.(fieldname2).sm_log=DATSTRUCT.(fieldname2).sm_log(inds,:);


%%

% Combine stim and stim catch into one matrix
sm=StatsStruct.(fieldname1).sm;
sm_log=log(sm);

DATSTRUCT.(fieldname1).sm_log= sm_log';

DATSTRUCT.All.sm_log= [DATSTRUCT.(fieldname2).sm_log; DATSTRUCT.(fieldname1).sm_log];


% total trials
NumAllTrials=size(DATSTRUCT.All.sm_log,1);

% PLOT
figure; hold on;
imagesc(Params.tf_bins.tSM,1:NumAllTrials,DATSTRUCT.All.sm_log)
colormap('spring');
title('amplitudes, sorted by stim onset');
xlabel('time');
ylabel('trial #')
colorbar

for i=1:length(TrigTimes_fromEpochStart);
    plot(TrigTimes_fromEpochStart(i)/1000,i,'ks');
end


%% RUNNING AVERAGE FOR AMPLITUDE (OVER TRIALS)    

% SMOOTH CATCH TRIALS
fieldname=fieldname1;

NumTimeBins=size(DATSTRUCT.(fieldname).sm_log,2);

for i=1:NumTimeBins; % for all time bins
    DATSTRUCT.(fieldname).trialsmoothed.sm_log(:,i)=smooth(DATSTRUCT.(fieldname).sm_log(:,i),Params.HEATMAP.runningbin); % smooth over trials
   
end

    % remove start and end as smooth binsize lowers (until is 1 at
    % endpoints)
trialstoremove=floor(Params.HEATMAP.runningbin/2); % how many from start and end
    
DATSTRUCT.(fieldname).trialsmoothed.sm_log(1:trialstoremove,:)=[];
DATSTRUCT.(fieldname).trialsmoothed.sm_log(end-trialstoremove+1:end,:)=[];

% % PLOT
% N=size(DATSTRUCT.(fieldname).trialsmoothed.sm_log,1); % num trials
% 
% % clims=[-3 3];
% figure; hold on;
% imagesc(Params.tf_bins.tSM,1:N,DATSTRUCT.(fieldname).trialsmoothed.sm_log)
% colormap('spring');
% title('Amplitude, sorted by stim onset and smoothed');
% xlabel('time');
% ylabel('trial #')
% colorbar


% FOR NOT CATCH
fieldname=fieldname2;

NumTimeBins=size(DATSTRUCT.(fieldname).sm_log,2);

for i=1:NumTimeBins; % for all time bins
    DATSTRUCT.(fieldname).trialsmoothed.sm_log(:,i)=smooth(DATSTRUCT.(fieldname).sm_log(:,i),Params.HEATMAP.runningbin); % smooth over trials  
end

    % remove start and end as smooth binsize lowers (until is 1 at
    % endpoints)
trialstoremove=floor(Params.HEATMAP.runningbin/2); % how many from start and end
    
DATSTRUCT.(fieldname).trialsmoothed.sm_log(1:trialstoremove,:)=[];
DATSTRUCT.(fieldname).trialsmoothed.sm_log(end-trialstoremove+1:end,:)=[];

% % PLOT
% 
% % clims=[-3 3];
% figure; hold on;
% imagesc(Params.tf_bins.tSM,1:N,DATSTRUCT.(fieldname).trialsmoothed.sm_log)
% colormap('spring');
% title('Stim trial Amplitude, sorted by stim onset and smoothed');
% xlabel('time');
% ylabel('trial #')
% colorbar
% 
% for i=1:N;
%     plot(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed(i)/1000,i,'kx');
% end

    
% COMBINE, and plot
% combine
DATSTRUCT.All.trialsmoothed.sm_log=[DATSTRUCT.(fieldname2).trialsmoothed.sm_log; DATSTRUCT.(fieldname1).trialsmoothed.sm_log];
N2=size(DATSTRUCT.All.trialsmoothed.sm_log,1); % num trials

% plot
figure; hold on;
imagesc(Params.tf_bins.tSM,1:N2,DATSTRUCT.All.trialsmoothed.sm_log)
colormap('spring');
title('Amplitude heat map, sorted by stim onset and smoothed');
xlabel('time');
ylabel('trial #')
colorbar

N=length(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed); % num trials
for i=1:N;
    plot(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed(i)/1000,i,'kx');
end

    
%% REALIGN DATA BY ANOTHER ONSET
% FIRST, use this to plot some to look at desired onset
figure; hold on;
plot(DATSTRUCT.(fieldname1).sm_log(37,:))
plot(DATSTRUCT.(fieldname2).sm_log(37,:))

figure; hold on;
plot(DATSTRUCT.(fieldname1).sm_log')
plot(DATSTRUCT.(fieldname2).sm_log')


% INPUTS
% where to find crossing?
if ~isfield(Params.HEATMAP,'risecross');
Params.HEATMAP.risecross=input('what amplitude to cross? '); % 1st cross rising amplitude
end

if ~isfield(Params.HEATMAP,'tmin');
Params.HEATMAP.tmin=input('start of window to look in? (in samples)'); % index at start of window of interest
end

if ~isfield(Params.HEATMAP,'tmax');
Params.HEATMAP.tmax=input('end of window?');
end


% RUN
list={fieldname1, fieldname2};

for j=1:length(list);
    fieldname=list{j};
    


% THIRD, get crossing values
NumTrials=size(DATSTRUCT.(fieldname).sm_log,1);
risecross=Params.HEATMAP.risecross;
tmin=Params.HEATMAP.tmin;
tmax=Params.HEATMAP.tmax;

for i=1:NumTrials;
    [~,ind]=min(abs(risecross-DATSTRUCT.(fieldname).sm_log(i,tmin:tmax)));
    Params.HEATMAP.NewAmplOnset_ind.(fieldname).onset(i)=ind+tmin; % new onset ind (relative to epoch start) in samples
    
    % convert to ms from start
    Params.HEATMAP.NewAmplOnset_ind.(fieldname).onset_ms(i)=Params.tf_bins.tSM(ind+tmin)*1000; % timebin
end
end


%% automatically get predur and postdur

EarliestOn=min([Params.HEATMAP.NewAmplOnset_ind.(fieldname1).onset Params.HEATMAP.NewAmplOnset_ind.(fieldname2).onset]);
LatestOn=max([Params.HEATMAP.NewAmplOnset_ind.(fieldname1).onset Params.HEATMAP.NewAmplOnset_ind.(fieldname2).onset]);

predur=EarliestOn-50; % in samples, 
songdur= size(DATSTRUCT.(fieldname).sm_log,2); % in samples
postdur=songdur-LatestOn-50;


Params.HEATMAP.NewAmplOnset_ind.StimCatch.predur_samps=predur;
Params.HEATMAP.NewAmplOnset_ind.StimCatch.postdur_samps=postdur;
Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.predur_samps=predur;
Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.postdur_samps=postdur;


%% PLOT AMPLITUDE, but with new onset.


list={fieldname1, fieldname2};

for j=1:length(list);
    fieldname=list{j};
    NumTrials=length(Params.HEATMAP.NewAmplOnset_ind.(fieldname).onset);
    
    
    % 1) compile all vectors
    SmLog_NewOnset=[];
    for i=1:NumTrials;
        onset=Params.HEATMAP.NewAmplOnset_ind.(fieldname).onset(i);
        t1=onset-Params.HEATMAP.NewAmplOnset_ind.(fieldname).predur_samps;
        t2=onset+Params.HEATMAP.NewAmplOnset_ind.(fieldname).postdur_samps;
        
        SmLog_NewOnset=[SmLog_NewOnset; DATSTRUCT.(fieldname).sm_log(i,t1:t2)];
        
    end
    DATSTRUCT.(fieldname).SmLog_NewOnset=SmLog_NewOnset;
end

% % 2) Convert stim times relative to new epoch onset
Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed= ...
    Params.HEATMAP.TrigTimes_fromEpochStart_delayed-(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.onset_ms-...
    Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.predur_samps*median(diff(Params.tf_bins.tSM))*1000); % i.e. is old trig time (ms) - new onset (ms) + new predur (ms)

% Get new tbins
NumSamps=size(DATSTRUCT.(fieldname2).SmLog_NewOnset,2);
binsize=median(diff(Params.tf_bins.tSM)); % in sec
Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM=0:binsize:NumSamps*binsize; % know number of samps. and binsize

    
%% COMBINE BOTH TRIAL TYPES
DATSTRUCT.All.SmLog_NewOnset= [DATSTRUCT.(fieldname2).SmLog_NewOnset; DATSTRUCT.(fieldname1).SmLog_NewOnset];

% total trials
NumAllTrials=size(DATSTRUCT.All.SmLog_NewOnset,1);


% PLOT
figure; hold on;
imagesc(Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM,1:NumAllTrials,DATSTRUCT.All.SmLog_NewOnset)
colormap('spring');
title('Amplitude, Locked to new syl onset');
xlabel('time');
ylabel('trial #')
colorbar


for i=1:length(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed);
    plot(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed(i)/1000,i,'ks');
end







%% RESORT BASED ON STIM ONSET
[Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted, inds]=sort(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed);
DATSTRUCT.(fieldname2).SmLog_NewOnset_sorted=DATSTRUCT.(fieldname2).SmLog_NewOnset(inds,:);


% PLOT
DATSTRUCT.All.SmLog_NewOnset_sorted= [DATSTRUCT.(fieldname2).SmLog_NewOnset_sorted; DATSTRUCT.(fieldname1).SmLog_NewOnset];


figure; hold on;
imagesc(Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM,1:NumAllTrials,DATSTRUCT.All.SmLog_NewOnset_sorted)
colormap('spring');
title('Amplitude, Locked to new syl onset - SORTED BY STIM');
xlabel('time');
ylabel('trial #')
colorbar


for i=1:NumStimTrials;
    plot(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted(i)/1000,i,'ks');
end



%% SMOOTH ACROSS TRIALS


% SMOOTH CATCH TRIALS
fieldname=fieldname1;

NumTimeBins=size(DATSTRUCT.(fieldname).SmLog_NewOnset,2);

for i=1:NumTimeBins; % for all time bins
    DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset(:,i)=smooth(DATSTRUCT.(fieldname).SmLog_NewOnset(:,i),Params.HEATMAP.runningbin); % smooth over trials
end

    % remove start and end as smooth binsize lowers (until is 1 at
    % endpoints)
trialstoremove=floor(Params.HEATMAP.runningbin/2); % how many from start and end
    
DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset(1:trialstoremove,:)=[];
DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset(end-trialstoremove+1:end,:)=[];



% FOR NOT CATCH
fieldname=fieldname2;

NumTimeBins=size(DATSTRUCT.(fieldname).SmLog_NewOnset_sorted,2);

for i=1:NumTimeBins; % for all time bins
    DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset_sorted(:,i)=smooth(DATSTRUCT.(fieldname).SmLog_NewOnset_sorted(:,i),Params.HEATMAP.runningbin); % smooth over trials  
end

    Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted_smoothed=...
        smooth(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted,Params.HEATMAP.runningbin);

    % remove start and end as smooth binsize lowers (until is 1 at
    % endpoints)
trialstoremove=floor(Params.HEATMAP.runningbin/2); % how many from start and end
    
DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset_sorted(1:trialstoremove,:)=[];
DATSTRUCT.(fieldname).trialsmoothed.SmLog_NewOnset_sorted(end-trialstoremove+1:end,:)=[];
Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted_smoothed(1:trialstoremove,:)=[];
Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted_smoothed(end-trialstoremove+1:end,:)=[];

    
% COMBINE, and plot
% combine
DATSTRUCT.All.trialsmoothed.SmLog_NewOnset_sorted=[DATSTRUCT.(fieldname2).trialsmoothed.SmLog_NewOnset_sorted; DATSTRUCT.(fieldname1).trialsmoothed.SmLog_NewOnset];
N2=size(DATSTRUCT.All.trialsmoothed.SmLog_NewOnset_sorted,1); % num trials

% plot
figure; hold on;
imagesc(Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM,1:N2,DATSTRUCT.All.trialsmoothed.SmLog_NewOnset_sorted)
colormap('spring');
title(['Amplitude heat map, smoothed, locked to onset of 1st syl binsize: ' Params.HEATMAP.runningbin]);
xlabel('time');
ylabel('trial #')
colorbar


N=length(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted_smoothed);
for i=1:N;
    plot(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed_sorted_smoothed(i)/1000,i,'ks');
end


%% Get mean amplitude, using data locked to onset
list={'StimCatch','StimNotCatch'};
PlotCols=lt_make_plot_colors(2,0);


    figure; hold on;
for i=1:length(list);
    fieldname=list{i};
    DATSTRUCT.(fieldname).SmLog_NewOnset_mean=mean(DATSTRUCT.(fieldname).SmLog_NewOnset,1);
    DATSTRUCT.(fieldname).SmLog_NewOnset_median=median(DATSTRUCT.(fieldname).SmLog_NewOnset,1);
    DATSTRUCT.(fieldname).SmLog_NewOnset_std=std(DATSTRUCT.(fieldname).SmLog_NewOnset,0,1);
    N=size(DATSTRUCT.(fieldname).SmLog_NewOnset,1);
    DATSTRUCT.(fieldname).SmLog_NewOnset_sem=DATSTRUCT.(fieldname).SmLog_NewOnset_std./sqrt(N-1);

    
    % PLOT
%     =plot(Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM(2:end),DATSTRUCT.(fieldname).SmLog_NewOnset_mean,'Color',PlotCols{i});
   
    shadedErrorBar(Params.HEATMAP.tf_bins_new.NewAmplOnset.tSM(2:end),DATSTRUCT.(fieldname).SmLog_NewOnset_mean,DATSTRUCT.(fieldname).SmLog_NewOnset_sem,...
        {'Linewidth',2,'Color',PlotCols{i}},1);
end
    





%% TRIGGER ON STIM

% automatically get pre and post durations
binsize=median(diff(Params.tf_bins.tSM)); % in sec
Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples=(1/1000)*Params.HEATMAP.TrigTimes_fromEpochStart_delayed./binsize;
Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples=round(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples);

EarliestOn=min(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples); % in samples
LatestOn=max(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples); % in samples

predur=EarliestOn-50; % in samples, 
songdur= size(DATSTRUCT.(fieldname).sm_log,2); % in samples
postdur=songdur-LatestOn-50;


Params.HEATMAP.TrigOnset.predur_samps=predur;
Params.HEATMAP.TrigOnset.postdur_samps=postdur;
Params.HEATMAP.TrigOnset.predur_ms=1000*binsize*predur;
Params.HEATMAP.TrigOnset.postdur_ms=1000*binsize*postdur;

% what is median trig onset, to use with non-trig trials
Params.HEATMAP.TrigOnset.MednTrigForStimCatch_samps=median(Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples);
Params.HEATMAP.TrigOnset.MednTrigForStimCatch_ms=median(Params.HEATMAP.TrigTimes_fromEpochStart_delayed);



% PLOT

% 1) not catch
NumTrials=length(Params.HEATMAP.NewAmplOnset_ind.(fieldname2).onset);
    
    
    % compile all vectors
    SmLog_StimOnset=[];
    for i=1:NumTrials;
        onset=Params.HEATMAP.TrigTimes_fromEpochStart_delayed_samples(i);
        t1=onset-Params.HEATMAP.TrigOnset.predur_samps;
        t2=onset+Params.HEATMAP.TrigOnset.postdur_samps;
        t1=round(t1);
        t2=round(t2);
        SmLog_StimOnset=[SmLog_StimOnset; DATSTRUCT.(fieldname2).sm_log(i,t1:t2)];
        
    end
    DATSTRUCT.(fieldname2).SmLog_StimOnset=SmLog_StimOnset;

% 2) catch (i.e. use mediate stim trig time)
NumTrials=length(Params.HEATMAP.NewAmplOnset_ind.(fieldname1).onset);

    % compile all vectors
    onset=Params.HEATMAP.TrigOnset.MednTrigForStimCatch_samps;
    predur=Params.HEATMAP.TrigOnset.predur_samps;
    postdur=Params.HEATMAP.TrigOnset.postdur_samps;
    
            t1=round(onset-predur);
        t2=round(onset+postdur);
        
    DATSTRUCT.(fieldname1).SmLog_StimOnset=DATSTRUCT.(fieldname1).sm_log(:,t1:t2);
    
    

% Get new tbins
NumSamps=size(DATSTRUCT.(fieldname2).SmLog_StimOnset,2);
binsize=median(diff(Params.tf_bins.tSM)); % in sec
Params.HEATMAP.tf_bins_new.StimOnset.tSM=0:binsize:NumSamps*binsize; % know number of samps. and binsize

%% PLOT, with locked to stim
DATSTRUCT.All.SmLog_StimOnset= [DATSTRUCT.(fieldname2).SmLog_StimOnset   ; DATSTRUCT.(fieldname1).SmLog_StimOnset];

% total trials
NumAllTrials=size(DATSTRUCT.All.SmLog_StimOnset,1);


% PLOT
figure; hold on;
imagesc(Params.HEATMAP.tf_bins_new.StimOnset.tSM,1:NumAllTrials,DATSTRUCT.All.SmLog_StimOnset)
colormap('spring');
title('Amplitude, Locked to trigger time (line) - median trigger used for catch trials');
xlabel('time');
ylabel('trial #')
colorbar

line([Params.HEATMAP.TrigOnset.predur_ms/1000 Params.HEATMAP.TrigOnset.predur_ms/1000], ylim)

% for i=1:NumStimTrials;
%     plot(Params.HEATMAP.NewAmplOnset_ind.StimNotCatch.TrigTimes_fromEpochStart_delayed(i)/1000,i,'ks');
% end



%% SAVE

tstamp=lt_get_timestamp(0);
Params.SAVEFOLDERS.HEATMAP=[Params.PLOT_TimeWindow_savedir '/HEATMAP_' tstamp];
mkdir(Params.SAVEFOLDERS.HEATMAP);
cd(Params.SAVEFOLDERS.HEATMAP);

save('DATSTRUCT','DATSTRUCT');
save('Params','Params');

savemultfigs

end