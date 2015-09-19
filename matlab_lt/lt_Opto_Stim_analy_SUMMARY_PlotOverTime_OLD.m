%% TO DO:
% 1) get stim epochs, and do analyses on those individually


clear all; close all
%% 

DATSTRUCT=struct('data',[]);
PARAMS=struct('global',[],'indiv',[]);

%% INPUTS
% 1) Data directories
BirdDir='/bluejay3/lucas/birds/wh73pk61/';

ListOfDirs=...
    {'031215_HVCChR2_XStimON_PitchShiftOFF_day4/lt_Opto_Stim_analy_batch.labeled.all_13Mar2015_1208X/PLOT_StimCatch_StimNotCatch_13Mar2015_1210/TimeWindow_13Mar2015_1213/OverTime_13Mar2015_1217',...
    '030815_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2009/PLOT_All_12Mar2015_2010/TimeWindow_12Mar2015_2010/OverTime_12Mar2015_2010',...
    '030815_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2010/PLOT_StimCatch_StimNotCatch_12Mar2015_2012/TimeWindow_12Mar2015_2012/OverTime_12Mar2015_2012',...
    '030915_HVCChR2_XStimOFF_PitchShiftON_day1/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2012/PLOT_All_12Mar2015_2014/TimeWindow_12Mar2015_2014/OverTime_12Mar2015_2014',...
    '031015_HVCChR2_XStimOFF_PitchShiftON_day2/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2014/PLOT_All_12Mar2015_2015/TimeWindow_12Mar2015_2015/OverTime_12Mar2015_2015',...
    '030715_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2015/PLOT_All_12Mar2015_2016/TimeWindow_12Mar2015_2016/OverTime_12Mar2015_2016',...
    '031115_HVCChR2_XStimOFF_PitchShiftON_day3/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2016/PLOT_All_12Mar2015_2016/TimeWindow_12Mar2015_2016/OverTime_12Mar2015_2016',...
    '030715_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2016/PLOT_StimCatch_StimNotCatch_12Mar2015_2018/TimeWindow_12Mar2015_2018/OverTime_12Mar2015_2018',...
    '031115_HVCChR2_XStimON_PitchShiftOFF_day3_actuallyoff/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2018/PLOT_StimCatch_StimNotCatch_12Mar2015_2021/TimeWindow_12Mar2015_2022/OverTime_12Mar2015_2022',...
    '030615_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2022/PLOT_StimCatch_StimNotCatch_12Mar2015_2025/TimeWindow_12Mar2015_2025/OverTime_12Mar2015_2025',...
    '031215_HVCChR2_XStimOFF_PitchShiftON_day4/lt_Opto_Stim_analy_batch.labeled.all_13Mar2015_1222X/PLOT_All_13Mar2015_1231/TimeWindow_13Mar2015_1232/OverTime_13Mar2015_1232',...
    '031415_HVCChR2_XStimON_PitchShiftOFF_day6/lt_Opto_Stim_analy_batch.labeled.all_15Mar2015_1139X/PLOT_StimCatch_StimNotCatch_15Mar2015_1143/TimeWindow_15Mar2015_1146/OverTime_15Mar2015_1146',...
    '031315_HVCChR2_XStimON_PitchShiftOFF_day5/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1248X/PLOT_StimCatch_StimNotCatch_14Mar2015_1252/TimeWindow_14Mar2015_1252/OverTime_14Mar2015_1253',...
    '031315_HVCChR2_XStimOFF_PitchShiftON_day5/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1511X/PLOT_All_14Mar2015_1512/TimeWindow_14Mar2015_1512/OverTime_14Mar2015_1513',...
    '031415_HVCChR2_XStimOFF_PitchShiftON_day6/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1517X/PLOT_All_14Mar2015_1518/TimeWindow_14Mar2015_1518/OverTime_14Mar2015_1518',...
    '031515_HVCChR2_XStimON_PitchShiftOFF_day7/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1215X/PLOT_StimCatch_StimNotCatch_17Mar2015_1217/TimeWindow_17Mar2015_1218/OverTime_17Mar2015_1218',...
    '031515_HVCChR2_XStimOFF_PitchShiftON_day7/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1057X/PLOT_All_17Mar2015_1059/TimeWindow_17Mar2015_1059/OverTime_17Mar2015_1100',...
    '031615_HVCChR2_XStimON_PitchShiftOFF_day8/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1112X/PLOT_StimCatch_StimNotCatch_17Mar2015_1115/TimeWindow_17Mar2015_1116/OverTime_17Mar2015_1117',...
    '031615_HVCChR2_XStimOFF_PitchShiftON_day8/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1128X/PLOT_All_17Mar2015_1130/TimeWindow_17Mar2015_1130/OverTime_17Mar2015_1130',...
    '031715_HVCChR2_XStimON_PitchShiftOFF_day9/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1137X/PLOT_StimCatch_StimNotCatch_17Mar2015_1138/TimeWindow_17Mar2015_1139/OverTime_17Mar2015_1140',...
    '031415_HVCChR2_XStimOFF_PitchShiftON_day6/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1152X/PLOT_All_17Mar2015_1154/TimeWindow_17Mar2015_1156/OverTime_17Mar2015_1156',...
    '031715_HVCChR2_XStimOFF_PitchShiftON_day9/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1623X/PLOT_All_17Mar2015_1624/TimeWindow_17Mar2015_1627/OverTime_17Mar2015_1627',...
    '031715_HVCChR2_XStimON_PitchShiftOFF_day9/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1024X/PLOT_StimCatch_StimNotCatch_18Mar2015_1027/TimeWindow_18Mar2015_1030/OverTime_18Mar2015_1032',...
    '031715_HVCChR2_XStimOFF_PitchShiftON_day9/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1039X/PLOT_All_18Mar2015_1040/TimeWindow_18Mar2015_1040/OverTime_18Mar2015_1040',...
    '031815_HVCChR2_XStimOFF_PitchShiftON_day10/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1439X/PLOT_All_18Mar2015_1440/TimeWindow_18Mar2015_1443/OverTime_18Mar2015_1443',...
    '031815_HVCChR2_XStimOFF_PitchShiftON_day10/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1257X/PLOT_All_19Mar2015_1259/TimeWindow_19Mar2015_1301/OverTime_19Mar2015_1301',...
    '031815_HVCChR2_XStimON_PitchShiftOFF_day10/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1244X/PLOT_StimCatch_StimNotCatch_19Mar2015_1245/TimeWindow_19Mar2015_1249/OverTime_19Mar2015_1250',...
    '031815_HVCChR2_DIRECTED/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1419X/PLOT_All_19Mar2015_1419/TimeWindow_19Mar2015_1439/OverTime_19Mar2015_1442'};
% these hold StatsStruct and Params, ideally in OverTime folder
% DONE:
% all to 3/11, from 3/6 (stim only)


% 2) Running average plots
PARAMS.global.RunBin=25; % number of renditions

% 3) What to plot
TimeFieldsOfInterest=[1 2 3 4];
statfield='ffvals'; % amplogvals, ffvals, entrvals

% 4) baseline days
PARAMS.global.BaselineDays=1:3;

%% Initialize
% put to PARAMS
PARAMS.global.BirdDir=BirdDir;
PARAMS.global.ListOfDirs=ListOfDirs;

%% Extract params

NumStructs=length(ListOfDirs);

%% First figure out range of dates
for i=1:NumStructs;
    
    dir=[BirdDir ListOfDirs{i}];
    
    prm=load([dir '/Params']);
    
    % get date
    slashes=strfind(prm.Params.savefolder,'/');
    dirname=prm.Params.savefolder(slashes(end-1)+1:slashes(end)-1);
    dateS{i}=dirname(1:6); % in mmddyy
end

date_nums=datenum(dateS,'mmddyy');

FirstDate_str=datestr(min(date_nums),'ddmmmyyyy');
FirstDate_num=min(date_nums);
LastDate_str=datestr(max(date_nums),'ddmmmyyyy');
LastDate_num=max(date_nums);


PARAMS.global.FirstDate_str=FirstDate_str;
PARAMS.global.FirstDate_num=FirstDate_num;
PARAMS.global.LastDate_str=LastDate_str;
PARAMS.global.LastDate_num=LastDate_num;

% Extract Structures and put into correct index based on date.
% FOR NOW, 0NLY DO PITCH - getting all values concatenated
% NOTE: only separates trials into 'All' (i.e. no stim) and 'StimCatch' and
% 'StimNotCatch' (with stim only ocurring when WN is off). puts all into
% one array per day.

c=1;
for i=1:NumStructs;
    
    dir=[BirdDir ListOfDirs{i}];
    
    sts=load([dir '/StatsStruct']);
    prm=load([dir '/Params']);
    
    % IF FIRST RUN, get some things for global params
    if c==1;
        PARAMS.global.TimeFields.char=prm.Params.TimeField;
        PARAMS.global.TimeFields.num=prm.Params.TimeWindowList;
        c=2;
    end
    
    % get date
    slashes=strfind(prm.Params.savefolder,'/');
    dirname=prm.Params.savefolder(slashes(end-1)+1:slashes(end)-1);
    date=dirname(1:6); % in ddmmyy
    dnum=datenum(date,'mmddyy');
    
    % FIGURE OUT IF THIS IS DIRECTED SONG EPOCH. if so, rename data as
    % such.
    if isfield(prm.Params,'DataInfo')==1;
        if strcmp(prm.Params.DataInfo,'Dir')==1;
           disp(['found dir song on ' date]);
    
            
        end
    end
    
    % Get info about struct
    % 1) index based on date
    ind=dnum-FirstDate_num+1;
    % 2) name of trial type field
    trialfields=fieldnames(sts.StatsStruct);
    % 3) how many fields?
    numtrialfields=length(trialfields);
    
    % Save Params
    try % if works, then day already has params, just add to it.
        tmp=PARAMS.indiv{ind};
        nparams=length(tmp);
        PARAMS.indiv{ind}.Params{nparams+1}=prm.Params;
    catch err % i.e. day is empty
        PARAMS.indiv{ind}.Params{1}=prm.Params;
    end
    
    
    % Save Data - slot into global struct
    for ii=1:length(prm.Params.TimeField); % for each timewindow
        twindfield=prm.Params.TimeField{ii};
        
        for iii=1:numtrialfields;
            
            tfield=trialfields{iii};
            
            % extract data
            ffvals=sts.StatsStruct.(tfield).WINDOWED.(twindfield).Pitch.vals;
            timevals=sts.StatsStruct.(tfield).WINDOWED.(twindfield).Time.hours;
            amplogvals=sts.StatsStruct.(tfield).WINDOWED.(twindfield).Ampl.vals_log;
            entrvals=sts.StatsStruct.(tfield).WINDOWED.(twindfield).WEntropy.vals;
            
            % slot in data:
            % CHECK IF DATA IS DIRECTED SONG:
            if isfield(prm.Params,'DataInfo')==1;
                if strcmp(prm.Params.DataInfo,'Dir')==1;
                    disp(['found dir song on ' date ' for trial type: ' tfield '; time window: ' twindfield]);
                    
                    tfield=['DIR_' tfield];
                end
            end
            
            
            try % if already exists, then add data to array.
                DATSTRUCT.data{ind}.(twindfield).(tfield).ffvals=[DATSTRUCT.data{ind}.(twindfield).(tfield).ffvals ffvals];
                DATSTRUCT.data{ind}.(twindfield).(tfield).timevals=[DATSTRUCT.data{ind}.(twindfield).(tfield).timevals timevals];
                DATSTRUCT.data{ind}.(twindfield).(tfield).amplogvals=[DATSTRUCT.data{ind}.(twindfield).(tfield).amplogvals amplogvals];
                DATSTRUCT.data{ind}.(twindfield).(tfield).entrvals=[DATSTRUCT.data{ind}.(twindfield).(tfield).entrvals entrvals];
                
                
            catch err
                % create new array
                DATSTRUCT.data{ind}.(twindfield).(tfield).ffvals=ffvals;
                DATSTRUCT.data{ind}.(twindfield).(tfield).timevals=timevals;
                DATSTRUCT.data{ind}.(twindfield).(tfield).amplogvals=amplogvals;
                DATSTRUCT.data{ind}.(twindfield).(tfield).entrvals=entrvals;
                
                
            end
            
        end
    end
end


NumDays=PARAMS.global.LastDate_num-PARAMS.global.FirstDate_num+1;




%% MAKE SURE ALL DATA are in temporal order (within each day)

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            
            % check if monotonically increasing timepoints
            timevals=DATSTRUCT.data{ii}.(timefield).(trialfield).timevals;
            
            if any(diff(timevals)<0)==1;
                % then they are not montonic up
                
                % sort
                [timevals,ind]=sort(timevals);
                statvals=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield)(ind); % this sorts
                
                DATSTRUCT.data{ii}.(timefield).(trialfield).timevals=timevals;
                DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield)=statvals;
            end
        end
    end
end


%% DETECT STIMULATION EPOCHS (assuming epochs are at least N hours apart)


%% COMPUTE RUNNING Median (each computed over individual epochs)

% for each day, get running average for each category (intracatergory)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            
            Yrun1=lt_running_stats(DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield),PARAMS.global.RunBin); % smooth data
            Yrun2=lt_running_stats(DATSTRUCT.data{ii}.(timefield).(trialfield).timevals,PARAMS.global.RunBin); % smooth time vals
            
            DATSTRUCT.data{ii}.(timefield).(trialfield).([statfield '_SmMedn'])=Yrun1.Median; % runing median
            DATSTRUCT.data{ii}.(timefield).(trialfield).([statfield '_SmMean'])=Yrun1.Mean; % runing median
            DATSTRUCT.data{ii}.(timefield).(trialfield).timevals_sm=Yrun2.Median;
        end
    end
end


%% PLOT smoothed
if (0)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};

    figure; hold on;
    title([timefield ' - ' statfield ' - smoothed with binsize: ' num2str(PARAMS.global.RunBin) ' rends.']);

    % plot each day
    for ii=1:NumDays;

        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);

        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            X=DATSTRUCT.data{ii}.(timefield).(trialfield).timevals_sm;
            X=X./24; % from hours to fraction of day.
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).([statfield '_SmMean']);

            switch trialfield
                case 'All'
                    plot(X+ii,Y,'k.','MarkerSize',6)
                case 'StimCatch'
                    plot(X+ii,Y,'r.','MarkerSize',6)
                case 'StimNotCatch'
                    plot(X+ii,Y,'g.','MarkerSize',6)
            end

%             % Plot mean and sem of this data type
%             ffmean=mean(Y);
%             ffstd=std(Y);
%             ffn=length(Y);
%             ffsem=ffstd/sqrt(ffn-1);
% 
%             switch trialfield
%                 case 'All'
%                     errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',9,'MarkerFaceColor','k')
%                 case 'StimCatch'
%                     errorbar(ii+1+iii*0.05,ffmean,ffsem,'ro','MarkerSize',9,'MarkerFaceColor','r')
%                 case 'StimNotCatch'
%                     errorbar(ii+1+iii*0.05,ffmean,ffsem,'go','MarkerSize',9,'MarkerFaceColor','g')
%             end
        end

%         % plot mean of the day (if there are multipel trial types)
%         if numtrialtypes>1;
%             Ytot=[];
%             for iii=1:numtrialtypes;
%                 trialfield=trialtypes{iii};
%                 Y=DATSTRUCT.data{ii}.(timefield).(trialfield).([statfield '_sm']);
%                 Ytot=[Ytot Y];
%             end
% 
%             ffmean=mean(Ytot);
%             ffstd=std(Ytot);
%             ffn=length(Ytot);
%             ffsem=ffstd/sqrt(ffn-1);
% 
%             errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
% 
%         end
    end
end

end
%% Compute running average for 2 classes: 1) stim, 2) stim catch + non-stim epoch.

% for each day, get running average for each category (intracatergory)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        % baseline vector
        if sum(strcmp(trialtypes,'All'))==1 && sum(strcmp(trialtypes,'StimCatch'))==1;
            % then combine baseline + stim catch into one vector
            Y_baseline=[DATSTRUCT.data{ii}.(timefield).All.(statfield) DATSTRUCT.data{ii}.(timefield).StimCatch.(statfield)];
            T_baseline=[DATSTRUCT.data{ii}.(timefield).All.timevals DATSTRUCT.data{ii}.(timefield).StimCatch.timevals];
            
        elseif sum(strcmp(trialtypes,'All'))==1 && sum(strcmp(trialtypes,'StimCatch'))==0;
            
            Y_baseline=DATSTRUCT.data{ii}.(timefield).All.(statfield);
            T_baseline=DATSTRUCT.data{ii}.(timefield).All.timevals;
            
        elseif sum(strcmp(trialtypes,'All'))==0 && sum(strcmp(trialtypes,'StimCatch'))==1;
            
            Y_baseline=DATSTRUCT.data{ii}.(timefield).StimCatch.(statfield);
            T_baseline=DATSTRUCT.data{ii}.(timefield).StimCatch.timevals;
        end
        
        % sort them so in temporal order
        [T_baseline, inds]=sort(T_baseline);
        Y_baseline=Y_baseline(inds);
        
        % put into DatStruct
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.(statfield)=Y_baseline;
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals=T_baseline;
        
        
        % SMOOTH
        Yrun1=lt_running_stats(DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.(statfield),PARAMS.global.RunBin); % smooth data
        Yrun2=lt_running_stats(DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals,PARAMS.global.RunBin); % smooth time vals
        
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.([statfield '_SmMedn'])=Yrun1.Median; % runing median
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.([statfield '_SmMean'])=Yrun1.Mean; % runing median
        
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals_sm=Yrun2.Median;
    end
end


%% PLOT - Raw and smoothed data.


% PLOT - all vals, along with running mean (combining stimcatch and
% baseline) (at median time points).
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ' - ' statfield ' - all renditions over days (running mean (time median))']);
    
    % plot each day
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            X=DATSTRUCT.data{ii}.(timefield).(trialfield).timevals;
            X=X./24; % from hours to fraction of day.
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield);
            
            switch trialfield
                case 'All'
                    plot(X+ii,Y,'k.')
                case 'StimCatch'
                    plot(X+ii,Y,'r.')
                case 'StimNotCatch'
                    plot(X+ii,Y,'g.')
            end
            
            % Plot mean and sem of this data type
            ffmean=mean(Y);
            ffstd=std(Y);
            ffn=length(Y);
            ffsem=ffstd/sqrt(ffn-1);
            
            switch trialfield
                case 'All'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',9,'MarkerFaceColor','k')
                case 'StimCatch'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'ro','MarkerSize',9,'MarkerFaceColor','r')
                case 'StimNotCatch'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'go','MarkerSize',9,'MarkerFaceColor','g')
            end
        end
        
        
        % PLOT running mean ((all and stim catch) vs. (stim not catch)
        % plot no-stim
        X=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals_sm;
        X=X./24; % from hours to fraction of day.
        Y=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.([statfield '_SmMean']);
        
        plot(X+ii,Y,'k-','LineWidth',2);
        
        % plot stim
        if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
            
            X=DATSTRUCT.data{ii}.(timefield).StimNotCatch.timevals_sm;
            X=X./24; % from hours to fraction of day.
            Y=DATSTRUCT.data{ii}.(timefield).StimNotCatch.([statfield '_SmMean']);
            
            plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
        end
        
        
        
        % plot mean of the day (if there are multipel trial types)
        if numtrialtypes>1;
            Ytot=[];
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield);
                Ytot=[Ytot Y];
            end
            
            ffmean=mean(Ytot);
            ffstd=std(Ytot);
            ffn=length(Ytot);
            ffsem=ffstd/sqrt(ffn-1);
            
            errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
            
        end
    end
end


% PLOT ONLY MEANS
NumDays=PARAMS.global.LastDate_num-PARAMS.global.FirstDate_num+1;
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ' - Mean ' statfield ' - all renditions over days']);
    
    % plot each day
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield);
            
            % Plot mean and sem of this data type
            ffmean=mean(Y);
            ffstd=std(Y);
            ffn=length(Y);
            ffsem=ffstd/sqrt(ffn-1);
            
            switch trialfield
                case 'All'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'ko','MarkerSize',9,'MarkerFaceColor','k')
                case 'StimCatch'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'ro','MarkerSize',9,'MarkerFaceColor','r')
                case 'StimNotCatch'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'go','MarkerSize',9,'MarkerFaceColor','g')
            end
        end
        
        % plot mean of the day (if there are multipel trial types)
        if numtrialtypes>1;
            Ytot=[];
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield);
                Ytot=[Ytot Y];
            end
            
            ffmean=mean(Ytot);
            ffstd=std(Ytot);
            ffn=length(Ytot);
            ffsem=ffstd/sqrt(ffn-1);
            
            errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
            
        end
    end
end


%% does effect of stim depend on recent slope of ff?
if (0)
% for each bin of stim epoch, take the last bin of non-stim and calculate
% 1) reversion and 2) slope of nonstim (this includes both baseline and
% nons-stim epochs

StimBinSizeList=[7 10 12 15];
% StimBinSizeList=[12]; % 12 seems good.

ii=8; % day
timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(4)};

for m=1:length(StimBinSizeList);
    
% perform with overlap.
StimBinSize=StimBinSizeList(m); % num renditions, (will get non-stim that are temporally closest to stim
% SlopeBinEdges=[-ceil((20-StimBinSize)/2) 20-ceil((20-StimBinSize)/2)]; % gives centered bins
SlopeBinEdges=[0 15]; % relative to 1st ind of stim bin (i.e. if [-3 7] then is 11 bins (i.e. [i-3:i+7], starting 3 inds before start of stim bin.
% SlopeBinEdges=[StimBinSize StimBinSize+20]; % relative to 1st ind of stim bin (i.e. if [-3 7] then is 11 bins (i.e. [i-3:i+7], starting 3 inds before start of stim bin.

MaxHoursApart=1; % don't want datasets with vals far apart.




% RUN - ----------------------------
% Get vectors
T_Stim=DATSTRUCT.data{ii}.(timefield).StimNotCatch.timevals;
Y_Stim=DATSTRUCT.data{ii}.(timefield).StimNotCatch.(statfield);

T_NoStim = DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals; % time vals
Y_NoStim = DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.(statfield); % stat

% Divvy up data based on stim vector and desired binsize (renditions)
% 1) how many bins will there be?
StimTrials=length(Y_Stim); % num stim trials
Nbins=StimTrials-StimBinSize+1;

% 2) for each bin, get info
c=1; % counter, since will throw out some bins
reversion=[];
slope=[];
medtime=[];
pval=[];
for k=1:Nbins;
    inds=k:k+StimBinSize-1;
    
    stimvals=Y_Stim(inds);
    stimtimes=T_Stim(inds);
    
    % skip analysis for bins with values that are too far apart (e.g. 1hr) - usually means different test epochs
    if max(abs(diff(stimtimes)))<MaxHoursApart;
        
        
        % nonstim vals, get times that are closest to stim - based on distance
        % from median timepoint of stim
        [~,nostimInds]=sort(abs(T_NoStim-median(stimtimes))); % gets inds of closest times to stim median
        
        indstokeep=nostimInds(1:StimBinSize); % the inds of the closest vals
        
        nostimvals=Y_NoStim(indstokeep); % get values for reversion
        nostimtimes=T_NoStim(indstokeep); % might not be in order.
        
        % sort nostimtimes
        [nostimtimes, ind]=sort(nostimtimes);
        nostimvals=nostimvals(ind);
        
        % get values for slope
        indstokeep_slope=min(indstokeep)+SlopeBinEdges(1):min(indstokeep)+SlopeBinEdges(2);
        
        % only continue if there is enough data to calculate slope
        if length(Y_NoStim)>=max(indstokeep_slope);
            
            nostimvals_slope=Y_NoStim(indstokeep_slope);
            nostimtimes_slope=T_NoStim(indstokeep_slope);
            
            % sort
            [nostimtimes_slope, ind]=sort(nostimtimes_slope);
            nostimvals_slope=nostimvals_slope(ind);
            
            
            
            % check to make sure nostim times are not too far apart.
            if max(abs(diff(nostimtimes)))<MaxHoursApart; % then continue
                
                % calculate reversion:
                reversion(c)=mean(stimvals)-mean(nostimvals);
                
                % get median stim time
                medtime(c)=median(stimtimes);
                
                % get slope (using no stim trials) using linear regression
                x=[ones(length(nostimtimes_slope),1) [1:length(nostimvals_slope)]'];
                y=nostimvals_slope';
                [b,bint,r,rint,stats]=regress(y, x);
                
                slope(c)=b(2);
                pval(c)=stats(3);
                
                % update counter
                c=c+1;
                
            end
        end
    end
end


% Correlation between reversion and slope
x=[ones(length(reversion),1) reversion'];
y=slope';
[b,bint,r,rint,stats]=regress(y,x);

% Plot 1 - 
figure; hold on;
title(['Slope of non-stim trials (binedges: ' num2str(SlopeBinEdges) ') vs. reversion caused by stim (binsize: ' num2str(StimBinSize) ')']);
xlabel('Stim minus Nonstim (FF, hz)');
ylabel('Slope of non-stim trials');

plot(reversion,slope,'k.');

% plot those with p<0.05 diff color
plot(reversion(pval<0.05),slope(pval<0.05),'r.');


Xlims=xlim;
Ylims=ylim;
X=linspace(Xlims(1),Xlims(2),length(reversion));
plot(X,b(1)+X.*b(2),'-');

htext=text(Xlims(1)+3,Ylims(1)+0.3,['R2 = ' num2str(stats(1)) '; p = ' num2str(stats(3))]);
set(htext,'FontSize',14,'Color','b')

% Plot 2, plot reversion and slope over time

end


% slope fluctuate up and down?
figure; hold on;
plot(medtime,(slope),'.b')
plot(medtime,(reversion),'.r')
line(xlim,[0 0]);

end
%% PLOT - OVERLAY ALL DAYS
if (0)
% overlay all days, smoothed plots (stim off vs. stim on).
NumDays=PARAMS.global.LastDate_num-PARAMS.global.FirstDate_num+1;
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on
    title([timefield ' - ' statfield ' - All days overlayed']);
    
    % plot each day
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).([statfield '_SmMean']);
            T=DATSTRUCT.data{ii}.(timefield).(trialfield).timevals_sm';

            
            switch trialfield
                case 'All'
                    plot(T,Y,'k.','MarkerSize',6)
                case 'StimCatch'
                    plot(T,Y,'r.','MarkerSize',6)
                case 'StimNotCatch'
                    plot(T,Y,'g.','MarkerSize',6)
            end
        end
    end
end
end


%% Bin all non-stim data by time (into hour bins)

HrBinEdges=0:1:24;

% TO ACCOUNT FOR DST - wake times being different from day to day
DST_mods{1}={'09Mar2015', '10Mar2015','11Mar2015','12Mar2015','15Mar2015'}; % format, celkl array, with col 1 being date, and col 2 being how much time to add to 7am to get actual wake time.
DST_mods{2}=[1,1,1,0.5,-1];


for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    % collect all data for each day
    Y=[]; % will fill with data for day;
    for ii=1:NumDays;
        Y=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.(statfield);
        T=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals;
        
        % Divide into 1hr bins.
        [~, bin]=histc(T,HrBinEdges); % for each value, tells you bin num (0hr is in bin 1, etc)
        
        % for all bins, put stat values into a cell, and get mean, etc
        for iii=1:length(HrBinEdges);
            HrBinnedData.data{ii}.(timefield).(statfield).ValsInHrBins{iii}=Y(bin==iii); % for each bin, put vals in
            HrBinnedData.data{ii}.(timefield).(statfield).Means(iii)=mean(Y(bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).STDs(iii)=std(Y(bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).N(iii)=sum((bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).SEM(iii)=std(Y(bin==iii))/sqrt(sum((bin==iii))-1);
        end
        HrBinnedData.data{ii}.(timefield).(statfield).HrBinEdges=HrBinEdges;
          
        
        % REDO, but adding or subtracting time for DST changes in lights on
        % time
        CurDate=datestr(FirstDate_num+ii-1,'ddmmmyyyy');
        
        if any(strcmp(DST_mods{1},CurDate)); % this date has diff wake time
            
            fudge=DST_mods{2}(strcmp(DST_mods{1},CurDate)); 
            T_fudge=T-fudge; % subtract since: if woke up at 8am (fudge = 1) then should move all data back by 1 hour
   
        % Divide into 1hr bins.
        [~, bin]=histc(T_fudge,HrBinEdges); % for each value, tells you bin num (0hr is in bin 1, etc)
        
        % for all bins, put stat values into a cell, and get mean, etc
        for iii=1:length(HrBinEdges);
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.ValsInHrBins{iii}=Y(bin==iii); % for each bin, put vals in
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.Means(iii)=mean(Y(bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.STDs(iii)=std(Y(bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.N(iii)=sum((bin==iii));
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.SEM(iii)=std(Y(bin==iii))/sqrt(sum((bin==iii))-1);
        end
            HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.T_fudge=T_fudge; % fudged time values.
        
        end
    end
end


%% WHAT IS BASELINE CIRCADIAN FLUCTUATION IN PITCH?

% 1) for all days, plot binned stats across time.  ask whether they show
% consistent pattern throughout day. if so, then subtract running avg from
% learning days above.

plotcols=lt_make_plot_colors(length(PARAMS.global.BaselineDays),0);
xplotmod=(HrBinEdges(2)-HrBinEdges(1))/2; % how much to add to x to get inbetween bin edges

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ', ' statfield ' - Baseline days, hr binned, with STD']);
    
    Ytot=[];
    X=HrBinnedData.data{ii}.(timefield).(statfield).HrBinEdges;
    
    for ii=PARAMS.global.BaselineDays;
        
        % for each baseline day, plot, overlayed, the binned data.
        Y=HrBinnedData.data{ii}.(timefield).(statfield).Means;
        Yerr=HrBinnedData.data{ii}.(timefield).(statfield).STDs;
        errorbar(X+rand*xplotmod,Y,Yerr,'o','Color',plotcols{ii},'MarkerFaceColor',plotcols{ii}); % x + 0.5 because want to plot between hours.
        
        % concat all baseline days to then take mean
        Ytot=[Ytot; Y];
    end
    
    xlim([6 23]);
    
    % take mean of all baseline days
    Ytot_mean=nanmean(Ytot,1);
    
    % PLOT mean and individual days
    figure; hold on;
    title([timefield ', ' statfield ', - Each bline day + mean across days']);
    plot(X+xplotmod,Ytot_mean,'k-o','MarkerFaceColor','k','MarkerSize',7);
    
    for ii=PARAMS.global.BaselineDays;
        plot(X+xplotmod,Ytot(ii,:),'b--o','MarkerFaceColor','b','MarkerSize',4);
    end
    xlim([6 23]);
    
    % OUTPUT
    % save baseline mean
    HrBinnedData.baseline.(timefield).(statfield).vals=Ytot;
    HrBinnedData.baseline.(timefield).(statfield).mean=Ytot_mean;
    HrBinnedData.baseline.(timefield).(statfield).HrBinEdges=X;
    
end

%% FOR ALL DATA (i.e. stim, stimcatch, all) subtract baseline from each val


for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    for ii=1:NumDays;
        
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            
            % data
            Yvals=DATSTRUCT.data{ii}.(timefield).(trialfield).(statfield);
            Tvals=DATSTRUCT.data{ii}.(timefield).(trialfield).timevals;
            
            % subtract baseline from data
            % 1) without DST fudge factor
            [~,bins]=histc(Tvals,HrBinEdges);
            
            Ymean_base=HrBinnedData.baseline.(timefield).(statfield).mean;
            tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
            
            Yvals_MinusBase=Yvals-tmp;
            
            % put back into structure
            DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.(statfield)=Yvals_MinusBase;
            
            
            % 2) with DST fudge (i.e. timevals).
            % first check if today has fudge. if so, fudge the Tvals:
            CurDate=datestr(FirstDate_num+ii-1,'ddmmmyyyy');
            
            if any(strcmp(DST_mods{1},CurDate)); % this date has diff wake time
                
                fudge=DST_mods{2}(strcmp(DST_mods{1},CurDate));
                T_fudge=Tvals-fudge; % subtract since: if woke up at 8am (fudge = 1) then should move all data back by 1 hour
                
                % using fudged Tvals:
                [~,bins]=histc(T_fudge,HrBinEdges);
                
                tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
                
                Yvals_MinusBase_fudge=Yvals-tmp;
                
                % put back into structure
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase_fudge;
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.timevals_fudgeDST=T_fudge;
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.actually_fudged='yes';
                
            else
                % fill it with non-fudged data
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase;
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.timevals_fudgeDST=Tvals;
                DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.actually_fudged='no';
                
            end
        end
        
        
        % ALSO DO FOR NON-stim (i.e. all + stim catch) trials
        % data
        Yvals=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.(statfield);
        Tvals=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals;
        
        % subtract baseline from data
        % 1) without DST fudge factor
        [~,bins]=histc(Tvals,HrBinEdges);
        
        Ymean_base=HrBinnedData.baseline.(timefield).(statfield).mean;
        tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
        
        Yvals_MinusBase=Yvals-tmp;
        
        % put back into structure
        DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.(statfield)=Yvals_MinusBase;
        
        
        % 2) with DST fudge (i.e. timevals).
        % first check if today has fudge. if so, fudge the Tvals:
        CurDate=datestr(FirstDate_num+ii-1,'ddmmmyyyy');
        
        if any(strcmp(DST_mods{1},CurDate)); % this date has diff wake time
            
            fudge=DST_mods{2}(strcmp(DST_mods{1},CurDate));
            T_fudge=Tvals-fudge; % subtract since: if woke up at 8am (fudge = 1) then should move all data back by 1 hour
            
            % using fudged Tvals:
            [~,bins]=histc(T_fudge,HrBinEdges);
            
            tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
            
            Yvals_MinusBase_fudge=Yvals-tmp;
            
            % put back into structure
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase_fudge;
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST=T_fudge;
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.actually_fudged='yes';
            
        else
            % fill it with non-fudged data
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase;
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST=Tvals;
            DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.actually_fudged='no';
        end
    end
end



% FOR ALL DAYS, PLOT RAW VALS
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ' - ' statfield ' - all renditions (baseline subtracted)']);
    
    % plot each day
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            X=DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.timevals_fudgeDST;
            X=X./24; % from hours to fraction of day.
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
            
            switch trialfield
                case 'All'
                    plot(X+ii,Y,'k.')
                case 'StimCatch'
                    plot(X+ii,Y,'r.')
                case 'StimNotCatch'
                    plot(X+ii,Y,'g.')
            end
            
            % Plot mean and sem of this data type
            ffmean=mean(Y);
            ffstd=std(Y);
            ffn=length(Y);
            ffsem=ffstd/sqrt(ffn-1);
            
            switch trialfield
                case 'All'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',9,'MarkerFaceColor','k')
                case 'StimCatch'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'ro','MarkerSize',9,'MarkerFaceColor','r')
                case 'StimNotCatch'
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'go','MarkerSize',9,'MarkerFaceColor','g')
            end
        end
        
        
        % PLOT running mean ((all and stim catch) vs. (stim not catch)
        % plot no-stim
%         X=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.timevals_sm;
%         X=X./24; % from hours to fraction of day.
%         Y=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.([statfield '_SmMean']);
%         
%         plot(X+ii,Y,'k-','LineWidth',2);
%         
%         % plot stim
%         if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
%             
%             X=DATSTRUCT.data{ii}.(timefield).StimNotCatch.timevals_sm;
%             X=X./24; % from hours to fraction of day.
%             Y=DATSTRUCT.data{ii}.(timefield).StimNotCatch.([statfield '_SmMean']);
%             
%             plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
%         end
%         
%         
        
        % plot mean of the day (if there are multipel trial types)
        if numtrialtypes>1;
            Ytot=[];
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                Ytot=[Ytot Y];
            end
            
            ffmean=mean(Ytot);
            ffstd=std(Ytot);
            ffn=length(Ytot);
            ffsem=ffstd/sqrt(ffn-1);
            
            errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
            
        end
    end
    line(xlim,[0 0]);
end


% PLOT ONLY MEANS
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ' - Mean ' statfield ' - all renditions over days']);
    
    % plot each day
    for ii=1:NumDays;
        
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            Y=DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
            
            % Plot mean and sem of this data type
            ffmean=mean(Y);
            ffstd=std(Y);
            ffn=length(Y);
            ffsem=ffstd/sqrt(ffn-1);
            
            switch trialfield
                case 'All'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'ko','MarkerSize',9,'MarkerFaceColor','k')
                case 'StimCatch'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'ro','MarkerSize',9,'MarkerFaceColor','r')
                case 'StimNotCatch'
                    errorbar(ii+iii*0.05,ffmean,ffsem,'go','MarkerSize',9,'MarkerFaceColor','g')
            end
        end
        
        % plot mean of the day (if there are multipel trial types)
        if numtrialtypes>1;
            Ytot=[];
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.(timefield).(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                Ytot=[Ytot Y];
            end
            
            ffmean=mean(Ytot);
            ffstd=std(Ytot);
            ffn=length(Ytot);
            ffsem=ffstd/sqrt(ffn-1);
            
            errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
            
        end
    end
        line(xlim,[0 0]);

end



%% DOES REVERSION CORRELATE WITH 1) amount of learning 2) slope preceding, 3) slope during stim epoch


% 1) LOOKING SIMPLY AT DAY MEANS
% first get amount of reversion per day

% collect learning and reversion
LearningTot=[];
LearningDiffTot=[];
ReversionTot=[];
ReversionForLearningDiff=[];
WithinDaySlopeTot=[];
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    hfig1=figure; hold on;
    title([timefield ' - learning (mean FF of nostim trials minus baseline (hr bins)) vs. reversion (stim - catch)']);
    xlabel('learning (hz)');
    ylabel('reversion (hz)');
    
    hfig2=figure; hold on;
    title([timefield ' - learning diff (today minus yesterday) vs. reversion (stim - catch)']);
    xlabel('learning Diff (hz)');
    ylabel('reversion (hz)');

    hfig3=figure; hold on;
    title([timefield ' - Within day slope (ffvals vs. time, nonimst) vs. reversion (stim - catch)']);
    xlabel('Non-stim slope (hz)');
    ylabel('reversion (hz)');
    
    
    
    % plot each day
    for ii=1:NumDays;
        
        if ~any(PARAMS.global.BaselineDays==ii); % ignore baseline days
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.(timefield));
            
            % only continue if today has both stim catch and notcatch
            if any(strcmp(trialtypes,'StimCatch')) && any(strcmp(trialtypes,'StimNotCatch'));
                
                % calculate reversion
                Ystim=mean(DATSTRUCT.data{ii}.(timefield).StimNotCatch.([statfield]));
                Ycatch=mean(DATSTRUCT.data{ii}.(timefield).StimCatch.([statfield]));
                
                REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.([statfield]) = Ystim-Ycatch;
                
                % collect all data
                ReversionTot=[ReversionTot REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.([statfield])];
                
                % PLOT
                % 1) Reversion vs. learning in that day
                learning = mean(DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])); % learning as mean of vals subtracted baseline
                
                figure(hfig1)
                plot(learning,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield),'ko','MarkerFaceColor','k');
                text(learning+1,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                
                % 2) change in learning (today minus yesterday);
                if ii~=1;
                    learningDiff=mean(DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']))...
                        -mean(DATSTRUCT_2.data{ii-1}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']));
                    LearningDiffTot=[LearningDiffTot learningDiff];
                    ReversionForLearningDiff=[ReversionForLearningDiff REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield)];
                    
                    figure(hfig2);
                    plot(learningDiff,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield),'ro','MarkerFaceColor','r');
                    text(learningDiff+1,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                end
                
                % 3) Within day slope (all non-stim trials).
                Tvals=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST;
                x=[ones(length(Tvals),1) Tvals']; % versus time
                y=DATSTRUCT_2.data{ii}.(timefield).All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])';
                [b,~,~,~,~]=regress(y,x); % get slope (b2)
                
                slope=b(2);
                
                figure(hfig3);
                plot(slope,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield),'go','MarkerFaceColor','g');
                text(slope+1,REVERSION.data{ii}.(timefield).DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                
                
                % collect all data
                LearningTot=[LearningTot learning];
                WithinDaySlopeTot=[WithinDaySlopeTot slope];
                
                
            end
        end
    end
    figure(hfig1);
    line(xlim,[0 0]);
    line([0 0],ylim);
    
    figure(hfig2);
    line(xlim,[0 0]);
    line([0 0],ylim);
   
end

% Plot all data combined
    %1) for learning
        x=[ones(length(LearningTot),1) LearningTot'];
y=ReversionTot';
[b,bint,r,rint,stats]=regress(y,x);

    figure; hold on;
    title('learning (mean FF of nostim trials minus baseline (hr bins)) vs. reversion (stim - catch) (All syls and days)');
    xlabel('learning (hz)');
    ylabel('reversion (hz)');
    plot(LearningTot,ReversionTot,'ko','MarkerFaceColor','k');
    line(xlim,[0 0]);
    line([0 0],ylim);
    Xlims=xlim;
    Ylims=ylim;
    
    plot(Xlims,b(1)+(Xlims).*b(2),'-k');
    htext=text(Xlims(1)+2,Ylims(1)+2,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))]);
    set(htext,'FontSize',14,'Color','b')

    
    % 2) for learning diff
    x=[ones(length(LearningDiffTot),1) LearningDiffTot'];
y=ReversionForLearningDiff';
[b,bint,r,rint,stats]=regress(y,x);

figure; hold on;
    title('Learning Diff vs. reversion (stim - catch) (All syls and days)');
    xlabel('learning diff (today - yesterday) (hz)');
    ylabel('reversion (hz)');
    plot(LearningDiffTot,ReversionForLearningDiff,'ro','MarkerFaceColor','r');
    line(xlim,[0 0]);
    line([0 0],ylim);
    Xlims=xlim;
    Ylims=ylim;
    
    plot(Xlims,b(1)+(Xlims).*b(2),'-k');
    htext=text(Xlims(1)+2,Ylims(1)+2,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))]);
    set(htext,'FontSize',14,'Color','b')

    
    % 3) within day slope
    x=[ones(length(WithinDaySlopeTot),1) WithinDaySlopeTot'];
y=ReversionTot';
[b,bint,r,rint,stats]=regress(y,x);

figure; hold on;
    title('Within day FF slope vs. reversion (stim - catch) (All syls and days)');
    xlabel('Within day slope (unit/hr)');
    ylabel('reversion (hz)');
    plot(WithinDaySlopeTot,ReversionTot,'go','MarkerFaceColor','g');
    line(xlim,[0 0]);
    line([0 0],ylim);
    Xlims=xlim;
    Ylims=ylim;
    
    plot(Xlims,b(1)+(Xlims).*b(2),'-k');
    htext=text(Xlims(1)+2,Ylims(1)+2,['R2=' num2str(stats(1)) '; p=' num2str(stats(3))]);
    set(htext,'FontSize',14,'Color','b')
    
    

%% TRANSFORM ALL hour binned data, by subtracting from baseline mean (for that time bin, mean of individual values).

if (0)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    Ybase=HrBinnedData.baseline.(timefield).(statfield).mean;
    
    for ii=1:NumDays;
        
        % determine if today is a DST modified wake time day.
        % if DST today, then use that timeshifted data.
        CurDate=datestr(FirstDate_num+ii-1,'ddmmmyyyy');
        
        if any(strcmp(DST_mods{1},CurDate)); % this date has diff wake time
            Y=HrBinnedData.data{ii}.(timefield).(statfield).FUDGED_DST.ValsInHrBins;
        else % today is not DST shifted
            Y=HrBinnedData.data{ii}.(timefield).(statfield).ValsInHrBins;
        end
        
        % for each hr bin, get values.
        for iii=1:length(Y); % hr bins
            
            % get vals
            HrBinnedData.data{ii}.(timefield).(statfield).MINUS_BASE_fudged.ValsInHrBins{iii}=Y{iii}-Ybase(iii); % each value in cell, subtract single number (i..e baseline mean of that bin)
            
            % get stats
            HrBinnedData.data{ii}.(timefield).(statfield).MINUS_BASE_fudged.Mean(iii)=mean(Y{iii}-Ybase(iii));
            HrBinnedData.data{ii}.(timefield).(statfield).MINUS_BASE_fudged.STDs(iii)=std(Y{iii}-Ybase(iii));
            HrBinnedData.data{ii}.(timefield).(statfield).MINUS_BASE_fudged.N(iii)=length(Y{iii}-Ybase(iii));
            HrBinnedData.data{ii}.(timefield).(statfield).MINUS_BASE_fudged.SEM(iii)=std(Y{iii}-Ybase(iii))/sqrt(length(Y{iii}-Ybase(iii))-1);
            
        end
    end
end



% PLOT ALL DAYS (without DST correction)
plotcols=lt_make_plot_colors(NumDays,0);
xplotmod=(HrBinEdges(2)-HrBinEdges(1))/2; % how much to add to x to get inbetween bin edges

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    title([timefield ', ' statfield ' - WN days, hr binned, with STD (not DST fudged)']);
    
    Ytot=[];
    X=HrBinnedData.data{ii}.(timefield).(statfield).HrBinEdges;
    
    for ii=PARAMS.global.BaselineDays(end)+1:NumDays;
        
        % for baseline day, plot, overlayed, the binned data.
        Y=HrBinnedData.data{ii}.(timefield).(statfield).Means;
        Yerr=HrBinnedData.data{ii}.(timefield).(statfield).STDs;
        errorbar(X+rand*xplotmod,Y,Yerr,'o','Color',plotcols{ii},'MarkerFaceColor',plotcols{ii}); % x + 0.5 because want to plot between hours.
        
        % concat all days to then take mean
        Ytot=[Ytot; Y];
    end
    
    xlim([6 23]);
    
    % take mean of all days
    Ytot_mean=nanmean(Ytot,1);
    
    % PLOT mean and individual days
    figure; hold on;
    title([timefield ', ' statfield ', - Each WN day + mean across days (mean of day means) (not DST fudged)']);
    plot(X+xplotmod,Ytot_mean,'k-o','MarkerFaceColor','k','MarkerSize',7);
    
    for ii=1:NumDays-PARAMS.global.BaselineDays(end);
        plot(X+xplotmod,Ytot(ii,:),'--o','Color',plotcols{ii},'MarkerFaceColor',plotcols{ii},'MarkerSize',4);
    end
    xlim([6 23]);
    
    % OUTPUT
%     % save global mean
%     HrBinnedData.baseline.(timefield).(statfield).vals=Ytot;
%     HrBinnedData.baseline.(timefield).(statfield).mean=Ytot_mean;
%     HrBinnedData.baseline.(timefield).(statfield).HrBinEdges=X;
    
end

end

%% FOR ALL STIM EPOCHS, plot stats:

% 1) Amount of reversion

% 2) Does reversion correlate with the recent slope of learning?









