function lt_seq_dep_pitch_Correlations(Params, AllDays_StructStatsStruct,DaysWanted);
%% LT 4/11/15 - Look at correlations between syllables
% TO DO:
% 1) How does correlation change thru learnig? - shoudl regress out
% learning. does that change relate to change in generalization?
% 2) Plot learning vs. correaltions (and structural similarity) - 
% 3) plot scatter of trials to get sense of strength of correlation.
% 4) Make sure corr and/or z-scoring is correct. got corr of 100 for using
% z-scores.
% 5) check that using song is a good proxy
% 6) plot sample sizes, etc.
% 7) seems to throw out a lot of songs. check that they are actually
% missing that data

% Inputs
% DaysWanted = 'baseline' % gets baseline corrs
% DaysWanted =10:13 % days 10:13

%% PARAMS

AllSylFields=fieldnames(AllDays_StructStatsStruct.IndivSyls);
BaselineDays=Params.SeqFilter.BaselineDays;


%% SUBTRACT DAY MEAN FOR EACH DAY

% all renditions, get z-score of feature vector relative to the
% mean/std feature vector of that day.
syllist=AllSylFields;
NumDays=datenum(Params.SeqFilter.LastDay)-datenum(Params.SeqFilter.FirstDay)+1;

for i=1:length(syllist);
    syl=syllist{i};
    
    % for each day, subtract mean of that day
    for ii=1:NumDays;
        
        inds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,ii);
        FVdayvals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(inds,:);
        FVdaymean=mean(FVdayvals);
        FVdaystd=std(FVdayvals);
        
        FVdayMinusMean=(FVdayvals-repmat(FVdaymean,size(FVdayvals,1),1));
        FVdayzscores=(FVdayvals-repmat(FVdaymean,size(FVdayvals,1),1))./repmat(FVdaystd,size(FVdayvals,1),1);
        
        
        AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(inds,:)=FVdayzscores;
        AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_MinusDayMean(inds,:)=FVdayMinusMean;
        
    end
end



%% PLOT ALL RENDITIONS TO VISUALIZE CORRELATIONS, ETC
for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder);

    syllist=Params.SeqFilter.SylLists.FieldsInOrder{j};

% 1) NOT SMOOTHED
figure; hold on;
title('All renditions - pitch vs. time');

PlotColors=lt_make_plot_colors(length(syllist),0);

for i=1:length(syllist);
    syl=syllist{i};
    
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot(i)=plot(times.FinalValue,vals,'-k');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot(i)=plot(times.FinalValue,vals,'.','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot(i)=plot(times.FinalValue,vals,'.','Color',rand*([1 0 0]));
      end
              
end

legend(hplot,syllist)

% 2) SMOOTHED
smbin=5; % rends to smooth over
figure; hold on;
title(['Smoothed (bin: ' num2str(smbin) ' rends): Pitch vs. time for all syls']);

for i=1:length(syllist);
    syl=syllist{i};
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    times_sm=lt_running_stats(times.FinalValue,smbin);
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
    vals_sm=lt_running_stats(vals,smbin);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-ok');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-o','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot(i)=plot(times_sm.Mean,vals_sm.Mean,'-o','Color',rand*([1 0 0]));
      end

end

legend(hplot,syllist)


% 3) Z-SCORED, same as above
% smoothed
smbin=5; % rends to smooth over
figure; hold on;
title(['Smoothed (bin: ' num2str(smbin) ' rends): Pitch vs. time for all syls']);
hfig=[];
for i=1:length(syllist);
    hfig(i)=subplot(length(syllist),1,i); hold on;
    syl=syllist{i};
    
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    times_sm=lt_running_stats(times.FinalValue,smbin);
    
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(:,2);
    vals_sm=lt_running_stats(vals,smbin);
   
      % PLOT
      % color based on relationship to targ syl.
      if  any(ismember(Params.SeqFilter.SylLists.TargetSyls,syl)); % target syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.k');
      elseif any(ismember(Params.SeqFilter.SylLists.SylsSame,syl)); % similar syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.','Color',rand*([0 1 0]));
      elseif any(ismember(Params.SeqFilter.SylLists.SylsDifferent,syl)); % similar syl
              hplot=plot(times_sm.Mean,vals_sm.Mean,'.','Color',rand*([1 0 0]));
      end

legend(hplot,syllist{i});
end

linkaxes(hfig,'x');

end


%% Get Correlation matrix by having song as a single value - 
% take advatange of fact that all renditions of all syllable times in a
% given song have exactly the same time values

% ORGANIZE ALL DATA INTO A MATRIX - each column is a syl, each row is a
% timepoint, and each element is average value in that timepoint (=song)

% convert to inds
if ischar(DaysWanted);
    if DaysWanted=='baseline';
        DayIndsWanted=BaselineDays;
    else
        disp('Problem, not sure what days wanted for corr matrix');
    end
elseif isnumeric(DaysWanted);
    DayIndsWanted=DaysWanted;
end



% get all times - these times required to contain all syls
syl=Params.SeqFilter.SylLists.TargetSyls{1}; % can use any syl to get times

% Only look at baseline days
bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DayIndsWanted); % inds of rends that are from baseline
times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(bInds); % times for all rends
times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times

SongTimes=unique(times.FinalValue); % find all unique times (i.e. songs), in baseline

% First concatenate syl names of all motifs (in order [motif 1 motif 2];
SylsAllMotifsOrdered={};
for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    SylsAllMotifsOrdered=[SylsAllMotifsOrdered Params.SeqFilter.SylLists.FieldsInOrder{j}];
end


% For each song time, get values for syls (1 mean value)
c=1;
syllist=SylsAllMotifsOrdered; % all syls from all motifs in order,

for i=1:length(SongTimes);
    curtime=SongTimes(i);
    
    for ii=1:length(syllist);
        syl=syllist{ii};
        
        % -- Extract times and values
        bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DayIndsWanted);
        times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(bInds); % times for all rends
        times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
        times=times.FinalValue;
        
        vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,2);
        
        % -- for each song time, get mean value
        % find values for that song
        inds=ismember(times,curtime);
        
        % Save info of how many rends in each song
        RendsPerSong(i,ii)=sum(inds);
        
        % if there are renditions, then get mean
        if any(inds)==0; % i.e. no value for this song
            MATRIX(c,ii)=nan;
        else
            MATRIX(c,ii)=mean(vals(inds));
        end
        
    end
    
    % check if any syls lacked data. if so, remove this row (by not
    % letting "c" add one more
    if any(isnan(MATRIX(c,:)));
        MATRIX(c,:)=[];
        disp(['THREW OUT SONG num (approx) ' num2str(c) ' because lacked data for at least one syl']);
        continue
    end
    
    c=c+1; % use this becuase some songs will throw out if does not have dat for all syls.
end

disp(['N (songs) used for corr: ' num2str(c-1)]);
disp(['N (songs) originally: ' num2str(length(SongTimes)) ]);


% GET CORRELATION MATRIX
[rho, p]=corr((MATRIX(:,:)));

% PLOT
figure; hold on;
colormap('jet');

subplot(2,1,1); hold on;
title(['Days ' num2str(DaysWanted) ': Correlation matrix (rho) betweens syls (mean val in song)']);
imagesc(rho,[-1 1]);
set(gca,'YTick',1:length(syllist));
set(gca,'YTickLabel',syllist);
set(gca,'XTick',1:length(syllist));
set(gca,'XTickLabel',syllist);

colorbar;

subplot(2,1,2); hold on;
title(['Days ' num2str(DaysWanted) ': matrix of p-values of correlation between syls (mean val in song)']);
imagesc(p,[-1 1]);
set(gca,'YTick',1:length(syllist));
set(gca,'YTickLabel',syllist);
set(gca,'XTick',1:length(syllist));
set(gca,'XTickLabel',syllist);

colorbar;




% TAKE 1D SLICE - correlations between target syl and others
for i=1:length(Params.SeqFilter.SylLists.TargetSyls); % for all targ syls
    targsyl=Params.SeqFilter.SylLists.TargetSyls{i};
    
    
    [rho, p]=corr(MATRIX); % gett correaltion matrix
    
    % find index of target syl.
    inds=strcmp(SylsAllMotifsOrdered,targsyl); % find where targ syl is in the corr matrix
    
    vals_rho=rho(:,inds); % extract column for targ syl vs. others
    vals_p=p(:,inds);
    
    figure; hold on;
    subplot(2,1,1); hold on;
    bar([vals_rho]);
    title(['Days ' num2str(DaysWanted) ' : rho of correlation between target (' targsyl ') and other syls']);
    set(gca,'XTick',1:length(SylsAllMotifsOrdered),'XTickLabel',SylsAllMotifsOrdered);
    
    subplot(2,1,2); hold on;
    bar([vals_p]);
    title(['Days ' num2str(DaysWanted) ' : p-value of corr between target (' targsyl ') and other syls']);
    set(gca,'XTick',1:length(SylsAllMotifsOrdered),'XTickLabel',SylsAllMotifsOrdered);
    
end

% -- CORRELATIONS PREDICT LEARNING?





%% CORRELATION FOR EACH DAY - look at evolution over time 
% OBSOLETE - USE ABOVE IN A FOR LOOP

%% GET CORRELATION MATRIX BY DIVIDING DATA USING TIME BINS
% will have a sliding time window - time windows in which both syls have
% enough samples will contribute to the correlation score.
HrBinSize=0.25; % size of bin, in hours
HrBinSlide=0.1; % amount to slide bin 
HrBinMinSamples=4; % minimum number of samples a bin must have (or it will be thrown out).







% Plot histogram of number of samples in all bins




