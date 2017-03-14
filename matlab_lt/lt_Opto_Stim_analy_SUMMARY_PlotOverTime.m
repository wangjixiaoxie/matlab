function lt_Opto_Stim_analy_SUMMARY_PlotOverTime(BirdDir,ListOfDirs,TimeFieldsOfInterest,statfield,BaselineDays,plotStimEpochs, MetaParams)
% plotStimEpochs=0 (skips plotting).
% MetaParams = data structure containing experiemnt info (e.g. what days
% paired stim + WN up).

%% TO DO:
% 1) get stim epochs, and do analyses on those individually

%% INPUTS examples
% 1) Data directories
% BirdDir='/bluejay3/lucas/birds/wh73pk61/';
%
% ListOfDirs=...
%     {'031215_HVCChR2_XStimON_PitchShiftOFF_day4/lt_Opto_Stim_analy_batch.labeled.all_13Mar2015_1208X/PLOT_StimCatch_StimNotCatch_13Mar2015_1210/TimeWindow_13Mar2015_1213/OverTime_13Mar2015_1217',...
%     '030815_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2009/PLOT_All_12Mar2015_2010/TimeWindow_12Mar2015_2010/OverTime_12Mar2015_2010',...
%     '030815_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2010/PLOT_StimCatch_StimNotCatch_12Mar2015_2012/TimeWindow_12Mar2015_2012/OverTime_12Mar2015_2012',...
%     '030915_HVCChR2_XStimOFF_PitchShiftON_day1/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2012/PLOT_All_12Mar2015_2014/TimeWindow_12Mar2015_2014/OverTime_12Mar2015_2014',...
%     '031015_HVCChR2_XStimOFF_PitchShiftON_day2/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2014/PLOT_All_12Mar2015_2015/TimeWindow_12Mar2015_2015/OverTime_12Mar2015_2015',...
%     '030715_HVCChR2_XStim_StimOFF/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2015/PLOT_All_12Mar2015_2016/TimeWindow_12Mar2015_2016/OverTime_12Mar2015_2016',...
%     '031115_HVCChR2_XStimOFF_PitchShiftON_day3/lt_Opto_Stim_analy_batch.labeled.catch_12Mar2015_2016/PLOT_All_12Mar2015_2016/TimeWindow_12Mar2015_2016/OverTime_12Mar2015_2016',...
%     '030715_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2016/PLOT_StimCatch_StimNotCatch_12Mar2015_2018/TimeWindow_12Mar2015_2018/OverTime_12Mar2015_2018',...
%     '031115_HVCChR2_XStimON_PitchShiftOFF_day3_actuallyoff/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2018/PLOT_StimCatch_StimNotCatch_12Mar2015_2021/TimeWindow_12Mar2015_2022/OverTime_12Mar2015_2022',...
%     '030615_HVCChR2_XStim_StimON_pulse/lt_Opto_Stim_analy_batch.labeled.all_12Mar2015_2022/PLOT_StimCatch_StimNotCatch_12Mar2015_2025/TimeWindow_12Mar2015_2025/OverTime_12Mar2015_2025',...
%     '031215_HVCChR2_XStimOFF_PitchShiftON_day4/lt_Opto_Stim_analy_batch.labeled.all_13Mar2015_1222X/PLOT_All_13Mar2015_1231/TimeWindow_13Mar2015_1232/OverTime_13Mar2015_1232',...
%     '031415_HVCChR2_XStimON_PitchShiftOFF_day6/lt_Opto_Stim_analy_batch.labeled.all_15Mar2015_1139X/PLOT_StimCatch_StimNotCatch_15Mar2015_1143/TimeWindow_15Mar2015_1146/OverTime_15Mar2015_1146',...
%     '031315_HVCChR2_XStimON_PitchShiftOFF_day5/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1248X/PLOT_StimCatch_StimNotCatch_14Mar2015_1252/TimeWindow_14Mar2015_1252/OverTime_14Mar2015_1253',...
%     '031315_HVCChR2_XStimOFF_PitchShiftON_day5/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1511X/PLOT_All_14Mar2015_1512/TimeWindow_14Mar2015_1512/OverTime_14Mar2015_1513',...
%     '031415_HVCChR2_XStimOFF_PitchShiftON_day6/lt_Opto_Stim_analy_batch.labeled.all_14Mar2015_1517X/PLOT_All_14Mar2015_1518/TimeWindow_14Mar2015_1518/OverTime_14Mar2015_1518',...
%     '031515_HVCChR2_XStimON_PitchShiftOFF_day7/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1215X/PLOT_StimCatch_StimNotCatch_17Mar2015_1217/TimeWindow_17Mar2015_1218/OverTime_17Mar2015_1218',...
%     '031515_HVCChR2_XStimOFF_PitchShiftON_day7/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1057X/PLOT_All_17Mar2015_1059/TimeWindow_17Mar2015_1059/OverTime_17Mar2015_1100',...
%     '031615_HVCChR2_XStimON_PitchShiftOFF_day8/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1112X/PLOT_StimCatch_StimNotCatch_17Mar2015_1115/TimeWindow_17Mar2015_1116/OverTime_17Mar2015_1117',...
%     '031615_HVCChR2_XStimOFF_PitchShiftON_day8/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1128X/PLOT_All_17Mar2015_1130/TimeWindow_17Mar2015_1130/OverTime_17Mar2015_1130',...
%     '031715_HVCChR2_XStimON_PitchShiftOFF_day9/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1137X/PLOT_StimCatch_StimNotCatch_17Mar2015_1138/TimeWindow_17Mar2015_1139/OverTime_17Mar2015_1140',...
%     '031415_HVCChR2_XStimOFF_PitchShiftON_day6/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1152X/PLOT_All_17Mar2015_1154/TimeWindow_17Mar2015_1156/OverTime_17Mar2015_1156',...
%     '031715_HVCChR2_XStimOFF_PitchShiftON_day9/lt_Opto_Stim_analy_batch.labeled.all_17Mar2015_1623X/PLOT_All_17Mar2015_1624/TimeWindow_17Mar2015_1627/OverTime_17Mar2015_1627',...
%     '031715_HVCChR2_XStimON_PitchShiftOFF_day9/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1024X/PLOT_StimCatch_StimNotCatch_18Mar2015_1027/TimeWindow_18Mar2015_1030/OverTime_18Mar2015_1032',...
%     '031715_HVCChR2_XStimOFF_PitchShiftON_day9/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1039X/PLOT_All_18Mar2015_1040/TimeWindow_18Mar2015_1040/OverTime_18Mar2015_1040',...
%     '031815_HVCChR2_XStimOFF_PitchShiftON_day10/lt_Opto_Stim_analy_batch.labeled.all_18Mar2015_1439X/PLOT_All_18Mar2015_1440/TimeWindow_18Mar2015_1443/OverTime_18Mar2015_1443',...
%     '031815_HVCChR2_XStimOFF_PitchShiftON_day10/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1257X/PLOT_All_19Mar2015_1259/TimeWindow_19Mar2015_1301/OverTime_19Mar2015_1301',...
%     '031815_HVCChR2_XStimON_PitchShiftOFF_day10/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1244X/PLOT_StimCatch_StimNotCatch_19Mar2015_1245/TimeWindow_19Mar2015_1249/OverTime_19Mar2015_1250',...
%     '031815_HVCChR2_DIRECTED/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1419X/PLOT_All_19Mar2015_1419/TimeWindow_19Mar2015_1439/OverTime_19Mar2015_1442',...
%     '031915_HVCChR2_DIRECTED/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1714X/PLOT_All_19Mar2015_1714/TimeWindow_19Mar2015_1715/OverTime_19Mar2015_1715',...
%     '031915_HVCChR2_XStimON_PitchShiftOFF_day11/lt_Opto_Stim_analy_batch.labeled.all_19Mar2015_1722X/PLOT_StimCatch_StimNotCatch_19Mar2015_1724/TimeWindow_19Mar2015_1739/OverTime_19Mar2015_1739',...
%     '032315_HVCChR2_NoteGroups_StimWnON_day3/lt_Opto_Stim_analy_batch.labeled.all_23Mar2015_1224X/PLOT_StimCatch_StimNotCatch_23Mar2015_1225/TimeWindow_23Mar2015_1225/OverTime_23Mar2015_1225',...
%     '032215_HVCChR2_NoteGroups_StimWnON_day2/lt_Opto_Stim_analy_batch.labeled.all_22Mar2015_1419X/PLOT_StimCatch_StimNotCatch_22Mar2015_1420/TimeWindow_22Mar2015_1421/OverTime_22Mar2015_1422',...
%     '032115_HVCChR2_NoteGroups_StimWnON_day1/lt_Opto_Stim_analy_batch.labeled.all_21Mar2015_2222X/PLOT_StimCatch_StimNotCatch_21Mar2015_2224/TimeWindow_21Mar2015_2226/OverTime_21Mar2015_2226'};
% these hold StatsStruct and Params, ideally in OverTime folder
% DONE:
% all to 3/11, from 3/6 (stim only)

%
% % 3) What to plot
% TimeFieldsOfInterest=[1 2 3 4];
% statfield='ffvals'; % amplogvals, ffvals, entrvals
%
% % 4) baseline days
% BaselineDays = 1:3
%


%% EXTRACT PARAMS

if ~exist('plotStimEpochs','var');
    plotStimEpochs=0;
end

DATSTRUCT=struct('data',[]);
PARAMS=struct('global',[],'indiv',[]);


PARAMS.global.BaselineDays=BaselineDays;
PARAMS.global.HrBinEdges=0:1:24; % bins to get baseline circadian flucturation.
PARAMS.global.RunBin=25; % number of renditions

% TO ACCOUNT FOR DST - wake times being different from day to day
PARAMS.global.DST_mods{1}={'09Mar2015', '10Mar2015','11Mar2015','12Mar2015','15Mar2015'}; % format, celkl array, with col 1 being date, and col 2 being how much time to add to 7am to get actual wake time.
PARAMS.global.DST_mods{2}=[1,1,1,0.5,-1];


%% Initialize
% put to PARAMS
PARAMS.global.BirdDir=BirdDir;
PARAMS.global.ListOfDirs=ListOfDirs;


NumStructs=length(ListOfDirs);

%% First figure out range of dates
for i=1:NumStructs;
    
    dir=[BirdDir ListOfDirs{i}];
    
    prm=load([dir '/Params']);
    
    % get date
    slashes=strfind(prm.Params.savefolder,'/');
    dirname=prm.Params.savefolder(slashes(5)+1:slashes(5)+6);
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
    dirname=prm.Params.savefolder(slashes(5)+1:slashes(5)+6);
    date=dirname(1:6); % in mmddyy
    dnum=datenum(date,'mmddyy');
    
    % Get info about struct
    % 1) index based on date
    ind=dnum-FirstDate_num+1;
    % 2) name of trial type field
    trialfields=fieldnames(sts.StatsStruct);
    
    
    
    % CHANGE NAMES TO FIT FOLLOWING CODE:
    % change NotStim to StimCatch
    if any(strcmp(trialfields,'NotStim'));
        % make sure StimCatch does not already exist;
        if ~any(strcmp(trialfields,'StimCatch'));
            % stimnotcatch is not there
            sts.StatsStruct.StimCatch=sts.StatsStruct.NotStim;
            sts.StatsStruct=rmfield(sts.StatsStruct,'NotStim');
            
        else
            % stimnotcatch already exists
            disp('problem - both NotStim and StimCatch in same day of structure - will overwrite, need to modify code to support this')
        end
    end
    
    % change Stim to StimNotCatch
    if any(strcmp(trialfields,'Stim'));
        % make sure StimNotCatch does not already exist;
        if ~any(strcmp(trialfields,'StimNotCatch'));
            % stimnotcatch is not there
            sts.StatsStruct.StimNotCatch=sts.StatsStruct.Stim;
            sts.StatsStruct=rmfield(sts.StatsStruct,'Stim');
        else
            % stimnotcatch already exists
            disp('problem - both Stim and StimNotCatch in same day of structure - will overwrite, need to modify code to support this')
        end
    end
    
    
    % Get current trial fields
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
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).ffvals=[DATSTRUCT.data{ind}.timewindow{ii}.(tfield).ffvals ffvals];
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).timevals=[DATSTRUCT.data{ind}.timewindow{ii}.(tfield).timevals timevals];
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).amplogvals=[DATSTRUCT.data{ind}.timewindow{ii}.(tfield).amplogvals amplogvals];
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).entrvals=[DATSTRUCT.data{ind}.timewindow{ii}.(tfield).entrvals entrvals];
                
                
            catch err
                % create new array
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).ffvals=ffvals;
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).timevals=timevals;
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).amplogvals=amplogvals;
                DATSTRUCT.data{ind}.timewindow{ii}.(tfield).entrvals=entrvals;
                
                
            end
            
        end
    end
end


NumDays=PARAMS.global.LastDate_num-PARAMS.global.FirstDate_num+1;


%% REMOVE OUTLIERS
% Empirically determined a value of 1 (for number of IQRs to use) as
% optimal, based on wh73pk61.
% only performs currently for the designated statfield (e.g. FF).
% if finds outlier in one time field, then removes that rendition from all
% time fields (makes compatible with stuff below that finds stim epochs
% assuming all time windows have same trial inds

for ii=1:NumDays;
    if ~isempty(DATSTRUCT.data{ii})
        % how many trial types today? e.g. StimCatch
        trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{1});
        numtrialtypes=length(trialtypes);
        
        for iii=1:numtrialtypes;
            trialfield=trialtypes{iii};
            
            Inds=[]; % collect inds across time fields
            for i=1:length(TimeFieldsOfInterest);
                
                Yvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                
                [Yvals_clean, outliers_hi outliers_lo]=lt_db_tukey_outlier_90tile(Yvals',1,1);
                
                % for plotting during troubleshooting
                %                 if ~isempty(outliers_hi)
                %                     figure; hold on; plot(Yvals,'o'), plot(outliers_hi, Yvals(outliers_hi),'or');
                %                 end
                %
                %                 if ~isempty(outliers_lo)
                %                     figure; hold on; plot(Yvals,'o'), plot(outliers_lo, Yvals(outliers_lo),'or');
                %                 end
                
                % what inds to remove?
                Inds=[Inds outliers_hi' outliers_lo'];
            end
            Inds=unique(Inds);
            
            % remove
            for i=1:length(TimeFieldsOfInterest);
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield)(Inds)=[];
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals(Inds)=[];
            end
        end
    end
end

% old version, ctreated each time window separately
% for i=1:length(TimeFieldsOfInterest); % for each window of interest
%     %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
%     for ii=1:NumDays;
%
%         if ~isempty(DATSTRUCT.data{ii})
%             % how many trial types today? e.g. StimCatch
%             trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
%             numtrialtypes=length(trialtypes);
%
%             for iii=1:numtrialtypes;
%                 trialfield=trialtypes{iii};
%
%                 Yvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
%
%                 [Yvals_clean, outliers_hi outliers_lo]=lt_db_tukey_outlier_90tile(Yvals',1,1);
%
%                 % for plotting during troubleshooting
%                 %                 if ~isempty(outliers_hi)
%                 %                     figure; hold on; plot(Yvals,'o'), plot(outliers_hi, Yvals(outliers_hi),'or');
%                 %                 end
%                 %
%                 %                 if ~isempty(outliers_lo)
%                 %                     figure; hold on; plot(Yvals,'o'), plot(outliers_lo, Yvals(outliers_lo),'or');
%                 %                 end
%
%                 % REPLACE data with outlier-removed data
%                 % what inds to remove?
%                 Inds=[outliers_hi' outliers_lo'];
%
%                 % remove
%                 DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield)(Inds)=[];
%                 DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals(Inds)=[];
%
%             end
%         end
%     end
% end



%% MAKE SURE ALL DATA are in temporal order (within each day)

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        
        if ~isempty(DATSTRUCT.data{ii})
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                
                % check if monotonically increasing timepoints
                timevals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals;
                
                if any(diff(timevals)<0)==1;
                    % then they are not montonic up
                    
                    % sort
                    [timevals,ind]=sort(timevals);
                    statvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield)(ind); % this sorts
                    
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals=timevals;
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield)=statvals;
                end
            end
        end
    end
end


%% COMBINE ALL NON-STIM DATA (i.e. stim catch and not-stim) AND SMOOTH that data

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            % combined no-stim vector
            if sum(strcmp(trialtypes,'All'))==1 && sum(strcmp(trialtypes,'StimCatch'))==1;
                % then combine baseline + stim catch into one vector
                Y_baseline=[DATSTRUCT.data{ii}.timewindow{i}.All.(statfield) DATSTRUCT.data{ii}.timewindow{i}.StimCatch.(statfield)];
                T_baseline=[DATSTRUCT.data{ii}.timewindow{i}.All.timevals DATSTRUCT.data{ii}.timewindow{i}.StimCatch.timevals];
                
            elseif sum(strcmp(trialtypes,'All'))==1 && sum(strcmp(trialtypes,'StimCatch'))==0;
                
                Y_baseline=DATSTRUCT.data{ii}.timewindow{i}.All.(statfield);
                T_baseline=DATSTRUCT.data{ii}.timewindow{i}.All.timevals;
                
            elseif sum(strcmp(trialtypes,'All'))==0 && sum(strcmp(trialtypes,'StimCatch'))==1;
                
                Y_baseline=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.(statfield);
                T_baseline=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.timevals;
            end
            
            % sort them so in temporal order
            [T_baseline, inds]=sort(T_baseline);
            Y_baseline=Y_baseline(inds);
            
            % put into DatStruct
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.(statfield)=Y_baseline;
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals=T_baseline;
        end
        
    end
end



%% DETECT STIMULATION EPOCHS (any stimulation that is preceded by a non-stim epoch)

% Tvals_diff_TOT=[]; % accumulates all diffs, across days
% for i=1:length(TimeFieldsOfInterest); % for each window of interest
%     %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
%     for ii=1:NumDays;
%         if ~isempty(DATSTRUCT.data{ii})
%
%             % how many trial types today? e.g. StimCatch
%             trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
%
%             if any(strcmp(trialtypes,'StimNotCatch')); % then this day has stim
%
%                 Tvals=[DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.timevals]; % combine stim catch and notcatch (i.e. assume they are in the same epoch).
%                 Tvals=sort(Tvals);
%
%                 Tvals_diff=diff(Tvals);
%
%                 Tvals_diff_TOT=[Tvals_diff_TOT Tvals_diff];
%
%             end
%
%         end
%     end
% end
%
%
% % plot histogram of Tval diffs
% edges=0:0.1:4;
% [n,~]=histc(Tvals_diff_TOT,edges);
%
% figure; hold on;
% bar(edges,log(n),'histc')

% STIM EPOCHS DEFINED AS ANY STIM PRECEDED BY NON-STIM.
% SKIP BASELINE DAYS
% just use the first time window (all time windows have same timings)
StimEpochs=cell(1, NumDays);
for ii=1:NumDays;
    
    if ~any(ii==BaselineDays); % skip baseline days
        
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            
            if any(strcmp(trialtypes,'StimNotCatch')); % then this day has stim
                
                if any(strcmp(trialtypes,'All')); % then this day also has non-stim
                    
                    %  make matrix with StimCatch times and id = 6
                    Tvals_StimCatch=[DATSTRUCT.data{ii}.timewindow{i}.StimCatch.timevals];
                    ID_StimCatch=6*ones(1,length(Tvals_StimCatch));
                    Index_StimCatch=1:length(Tvals_StimCatch);
                    StimCatchMat=[Index_StimCatch; ID_StimCatch; Tvals_StimCatch];
                    
                    %  make matrix with StimNotCatch times and id=5;
                    Tvals_StimNotCatch = [DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.timevals];
                    ID_StimNotCatch=5*ones(1,length(Tvals_StimNotCatch));
                    Index_StimNotCatch=1:length(Tvals_StimNotCatch);
                    StimNotCatchMat=[Index_StimNotCatch; ID_StimNotCatch; Tvals_StimNotCatch];
                    
                    % make matrix with nonstim times, and id = 0;
                    Tvals_notstim=DATSTRUCT.data{ii}.timewindow{i}.All.timevals;
                    ID_notstim=zeros(1,length(Tvals_notstim));
                    Index_NotStim=1:length(Tvals_notstim);
                    NotStimMat=[Index_NotStim; ID_notstim; Tvals_notstim];
                    
                    % combine and sort all combined time values
                    CombinedMat=[StimCatchMat StimNotCatchMat NotStimMat];
                    % sort
                    [~, ind]=sort(CombinedMat(3,:));
                    CombinedMat_sorted=CombinedMat(:,ind);
                    
                    % this gives history of trial types across time
                    % (transitions).
                    ID_diffs=diff(CombinedMat_sorted(2,:));
                    
                    % FIND TRANSITIONS
                    Nstim_to_Stim=find(ID_diffs>4);
                    Stim_to_Nstim=find(ID_diffs<-4);
                    
                    %                 NumStimEpochs=max([length(Nstim_to_Stim) length(Stim_to_Nstim)]);
                    
                    % extract indices of each stim epoch
                    
                    % 1) case 1
                    if length(Nstim_to_Stim)>length(Stim_to_Nstim);
                        for k=1:length(Nstim_to_Stim);
                            
                            % stim trials
                            if k==length(Nstim_to_Stim); % i.e. this epoch goes to end of day
                                EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k)+1:end);
                            else
                                EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k)+1:Stim_to_Nstim(k));
                            end
                            
                            StimEpochs{ii}.epoch{k}.StimNotCatch=EpochMat(1,EpochMat(2,:)==5); % Indices (within trial type) of trials that fall within this stim epoch
                            StimEpochs{ii}.epoch{k}.StimCatch=EpochMat(1,EpochMat(2,:)==6);
                            
                            
                            % notstim (all) trials preceding
                            if k==1;
                                EpochMat=CombinedMat_sorted(:,1:Nstim_to_Stim(k));
                            else
                                EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k-1)+1:Nstim_to_Stim(k));
                            end
                            StimEpochs{ii}.epoch{k}.All_preceding=EpochMat(1,EpochMat(2,:)==0);
                            
                            % notstim (all) trials post
                            if k~=length(Nstim_to_Stim);
                                EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:Nstim_to_Stim(k+1));
                                
                                StimEpochs{ii}.epoch{k}.All_post=EpochMat(1,EpochMat(2,:)==0);
                            end
                        end
                        
                        
                        
                        
                        
                        % 2) case 2
                    elseif length(Stim_to_Nstim)>length(Nstim_to_Stim);
                        
                        for k=1:length(Stim_to_Nstim);
                            % stim trials
                            if k==1;
                                EpochMat=CombinedMat_sorted(:,1:Stim_to_Nstim(k));
                            else
                                EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k-1)+1:Stim_to_Nstim(k));
                            end
                            
                            StimEpochs{ii}.epoch{k}.StimNotCatch=EpochMat(1,EpochMat(2,:)==5); % Indices (within trial type) of trials that fall within this stim epoch
                            StimEpochs{ii}.epoch{k}.StimCatch=EpochMat(1,EpochMat(2,:)==6);
                            
                            
                            % notstim (all) trials preceding
                            if k~=1;
                                EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k-1)+1:Nstim_to_Stim(k));
                                StimEpochs{ii}.epoch{k}.All_preceding=EpochMat(1,EpochMat(2,:)==0);
                            end
                            
                            
                            % notstim trials post
                            if k==length(Stim_to_Nstim);
                                EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:end);
                            else
                                EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:Nstim_to_Stim(k));
                            end
                            
                            StimEpochs{ii}.epoch{k}.All_post=EpochMat(1,EpochMat(2,:)==0);
                        end
                        
                        
                        % 3) case 3
                    elseif length(Stim_to_Nstim)==length(Nstim_to_Stim);
                        
                        % 3a) case 3a
                        if Stim_to_Nstim(1)<Nstim_to_Stim(1); % day starts with stim
                            
                            for k=1:length(Stim_to_Nstim)+1;
                                
                                % stim trials
                                if k==1;
                                    EpochMat=CombinedMat_sorted(:,1:Stim_to_Nstim(k));
                                else
                                    if k==length(Stim_to_Nstim)+1; % i.e. at the end
                                        EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k-1)+1:end);
                                    else
                                        EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k-1)+1:Stim_to_Nstim(k));
                                    end
                                end
                                StimEpochs{ii}.epoch{k}.StimNotCatch=EpochMat(1,EpochMat(2,:)==5); % Indices (within trial type) of trials that fall within this stim epoch
                                StimEpochs{ii}.epoch{k}.StimCatch=EpochMat(1,EpochMat(2,:)==6);
                                
                                % notstim (pre)
                                if k~=1
                                    EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k-1)+1:Nstim_to_Stim(k-1));
                                    StimEpochs{ii}.epoch{k}.All_preceding=EpochMat(1,EpochMat(2,:)==0);
                                end
                                
                                
                                % nostim (post)
                                if k~=length(Stim_to_Nstim)+1;
                                    EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:Nstim_to_Stim(k));
                                    StimEpochs{ii}.epoch{k}.All_post=EpochMat(1,EpochMat(2,:)==0);
                                    
                                end
                                
                            end
                            
                            % 3b) case 3b
                        else % then day starts with notstim
                            
                            for k=1:length(Stim_to_Nstim)
                                
                                % stim trials
                                EpochMat=CombinedMat_sorted(:,Nstim_to_Stim(k)+1:Stim_to_Nstim(k));
                                StimEpochs{ii}.epoch{k}.StimNotCatch=EpochMat(1,EpochMat(2,:)==5); % Indices (within trial type) of trials that fall within this stim epoch
                                StimEpochs{ii}.epoch{k}.StimCatch=EpochMat(1,EpochMat(2,:)==6);
                                
                                % notstim (pre)
                                if k==1
                                    EpochMat=CombinedMat_sorted(:,1:Nstim_to_Stim(k));
                                else
                                    EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k-1)+1:Nstim_to_Stim(k));
                                end
                                StimEpochs{ii}.epoch{k}.All_preceding=EpochMat(1,EpochMat(2,:)==0);
                                
                                
                                % nostim (post)
                                if k==length(Nstim_to_Stim);
                                    EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:end);
                                else
                                    EpochMat=CombinedMat_sorted(:,Stim_to_Nstim(k)+1:Nstim_to_Stim(k+1));
                                end
                                StimEpochs{ii}.epoch{k}.All_post=EpochMat(1,EpochMat(2,:)==0);
                                
                            end
                            
                        end
                    end
                    
                    
                    
                else % this day does not have non-stim, so all stim trials defined as contained within one stim epoch
                    % just a stim epoch, no pre and post
                    k=1;
                    StimEpochs{ii}.epoch{k}.StimNotCatch=1:length(DATSTRUCT.data{ii}.timewindow{1}.StimNotCatch.timevals);
                    StimEpochs{ii}.epoch{k}.StimCatch=1:length(DATSTRUCT.data{ii}.timewindow{1}.StimCatch.timevals);
                end
                
            end
            
        end
    end
end




%% SUBTRACT BASELINE FROM DATA (baseline determined in time bins to account for circadian)
% Bin all non-stim data by time (into hour bins)
HrBinEdges=PARAMS.global.HrBinEdges;

% TO ACCOUNT FOR DST - wake times being different from day to day
DST_mods=PARAMS.global.DST_mods;


for i=1:length(TimeFieldsOfInterest); % for each window of interest
    %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    % collect all data for each day
    Y=[]; % will fill with data for day;
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT_2.data{ii})
            
            Y=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.(statfield);
            T=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals;
            
            % Divide into bins.
            [~, bin]=histc(T,HrBinEdges); % for each value, tells you bin num (0hr is in bin 1, etc)
            
            % for all bins, put stat values into a cell, and get mean, etc
            for iii=1:length(HrBinEdges);
                HrBinnedData.data{ii}.timewindow{i}.(statfield).ValsInHrBins{iii}=Y(bin==iii); % for each bin, put vals in
                HrBinnedData.data{ii}.timewindow{i}.(statfield).Means(iii)=mean(Y(bin==iii));
                HrBinnedData.data{ii}.timewindow{i}.(statfield).STDs(iii)=std(Y(bin==iii));
                HrBinnedData.data{ii}.timewindow{i}.(statfield).N(iii)=sum((bin==iii));
                HrBinnedData.data{ii}.timewindow{i}.(statfield).SEM(iii)=std(Y(bin==iii))/sqrt(sum((bin==iii))-1);
            end
            HrBinnedData.data{ii}.timewindow{i}.(statfield).HrBinEdges=HrBinEdges;
            
            
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
                    HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.ValsInHrBins{iii}=Y(bin==iii); % for each bin, put vals in
                    HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.Means(iii)=mean(Y(bin==iii));
                    HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.STDs(iii)=std(Y(bin==iii));
                    HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.N(iii)=sum((bin==iii));
                    HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.SEM(iii)=std(Y(bin==iii))/sqrt(sum((bin==iii))-1);
                end
                HrBinnedData.data{ii}.timewindow{i}.(statfield).FUDGED_DST.T_fudge=T_fudge; % fudged time values.
                
            end
        end
    end
end


%% WHAT IS BASELINE CIRCADIAN FLUCTUATION IN PITCH?
% 1) for baseline days, plot binned stats across time.  ask whether they show
% consistent pattern throughout day. subtract those values from WN day vales

plotcols=lt_make_plot_colors(length(PARAMS.global.BaselineDays),0);
xplotmod=(HrBinEdges(2)-HrBinEdges(1))/2; % how much to add to x to get in between bin edges

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    figure; hold on;
    subplot(2,1,1);
    hold on; title([timefield ', ' statfield ' - Baseline days, hr binned, with STD']);
    
    Ytot=[];
    X=HrBinnedData.data{ii}.timewindow{i}.(statfield).HrBinEdges;
    
    for ii=PARAMS.global.BaselineDays;
        
        % for each baseline day, plot, overlayed, the binned data.
        Y=HrBinnedData.data{ii}.timewindow{i}.(statfield).Means;
        Yerr=HrBinnedData.data{ii}.timewindow{i}.(statfield).STDs;
        errorbar(X+rand*xplotmod,Y,Yerr,'o','Color',plotcols{ii},'MarkerFaceColor',plotcols{ii}); % x + 0.5 because want to plot between hours.
        
        % concat all baseline days to then take mean
        Ytot=[Ytot; Y];
    end
    
    xlim([6 23]);
    
    % take mean of all baseline days
    Ytot_mean=nanmean(Ytot,1);
    
    % PLOT mean and individual days
    subplot(2,1,2); hold on;
    title([timefield ', ' statfield ', - Each bline day + mean across days']);
    plot(X+xplotmod,Ytot_mean,'k-o','MarkerFaceColor','k','MarkerSize',7);
    
    for ii=PARAMS.global.BaselineDays;
        plot(X+xplotmod,Ytot(ii,:),'b--o','MarkerFaceColor','b','MarkerSize',4);
    end
    xlim([6 23]);
    
    % OUTPUT
    % save baseline mean
    HrBinnedData.baseline.timewindow{i}.(statfield).vals=Ytot;
    HrBinnedData.baseline.timewindow{i}.(statfield).mean=Ytot_mean;
    HrBinnedData.baseline.timewindow{i}.(statfield).HrBinEdges=X;
end

%% FOR ALL DATA (i.e. stim, stimcatch, all) subtract baseline from each val

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    for ii=1:NumDays;
        
        if ~isempty(DATSTRUCT.data{ii})
            
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                
                % data
                Yvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                Tvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals;
                
                % subtract baseline from data
                % 1) without DST fudge factor
                [~,bins]=histc(Tvals,HrBinEdges);
                
                Ymean_base=HrBinnedData.baseline.timewindow{i}.(statfield).mean;
                tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
                
                Yvals_MinusBase=Yvals-tmp;
                
                % put back into structure
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.(statfield)=Yvals_MinusBase;
                
                
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
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase_fudge;
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST=T_fudge;
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.actually_fudged='yes';
                    
                else
                    % fill it with non-fudged data
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase;
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST=Tvals;
                    DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.actually_fudged='no';
                    
                end
            end
            
            
            % ALSO DO FOR NON-stim (i.e. all + stim catch) trials
            % data
            Yvals=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.(statfield);
            Tvals=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals;
            
            % subtract baseline from data
            % 1) without DST fudge factor
            [~,bins]=histc(Tvals,HrBinEdges);
            
            Ymean_base=HrBinnedData.baseline.timewindow{i}.(statfield).mean;
            tmp=Ymean_base(bins); % converted to matrix with baseline values to subtract (length of data).
            
            Yvals_MinusBase=Yvals-tmp;
            
            % put back into structure
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.(statfield)=Yvals_MinusBase;
            
            
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
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase_fudge;
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST=T_fudge;
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.actually_fudged='yes';
                
            else
                % fill it with non-fudged data
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])=Yvals_MinusBase;
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST=Tvals;
                DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.actually_fudged='no';
            end
        end
    end
end


%% COMPUTE RUNNING avg (each computed over individual epochs)

% for each day, get running average for each category (intracatergory)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    %     timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    for ii=1:NumDays;
        
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                
                Yrun1=lt_running_stats(DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield),PARAMS.global.RunBin); % smooth data
                Yrun2=lt_running_stats(DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals,PARAMS.global.RunBin); % smooth time vals
                
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_SmMedn'])=Yrun1.Median; % runing median
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_SmMean'])=Yrun1.Mean; % runing median
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_SmSTD'])=Yrun1.STD;
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals_sm=Yrun2.Median;
                
                
                % DO SAME FOR BASELINE SUBTRACTED DATA
                Yrun1=lt_running_stats(DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']),PARAMS.global.RunBin); % smooth data
                Yrun2=lt_running_stats(DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST,PARAMS.global.RunBin); % smooth time vals
                
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST_SmMedn'])=Yrun1.Median;
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST_SmMean'])=Yrun1.Mean;
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD'])=Yrun1.STD;
                DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST_sm=Yrun2.Median;
                
            end
            
            % DO IT FOR NON-stim combined data
            % non-baseline subtracte data:
            Yrun1=lt_running_stats(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.(statfield),PARAMS.global.RunBin); % smooth data
            Yrun2=lt_running_stats(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals,PARAMS.global.RunBin); % smooth time vals
            
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.([statfield '_SmMedn'])=Yrun1.Median; % runing median
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.([statfield '_SmMean'])=Yrun1.Mean; % runing median
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.([statfield '_SmSTD'])=Yrun1.STD; % runing median
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals_sm=Yrun2.Median; % runing median
            
            
            % baseline subtractged data
            Yrun1=lt_running_stats(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']),PARAMS.global.RunBin);
            Yrun2=lt_running_stats(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST,PARAMS.global.RunBin); % smooth time vals
            
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST_SmMedn'])=Yrun1.Median;
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean'])=Yrun1.Mean;
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD'])=Yrun1.STD;
            DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST_sm=Yrun2.Median;
            
        end
    end
end


%% PLOT smoothed alone - IGNORE
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


%% PLOT - Raw and smoothed data. combined

if (0)
    % PLOT - all vals, along with running mean (combining stimcatch and
    % baseline) (at median time points).
    for i=1:length(TimeFieldsOfInterest); % for each window of interest
        timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
        
        figure; hold on;
        title([timefield ' - ' statfield ' - all renditions over days (running mean (time median))']);
        
        % plot each day
        for ii=1:NumDays;
            
            if ~isempty(DATSTRUCT.data{ii})
                
                % how many trial types today? e.g. StimCatch
                trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
                numtrialtypes=length(trialtypes);
                
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals;
                    X=X./24; % from hours to fraction of day.
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                    
                    switch trialfield
                        case 'All'
                            plot(X+ii,Y,'k.')
                        case 'StimCatch'
                            plot(X+ii,Y,'r.')
                        case 'StimNotCatch'
                            plot(X+ii,Y,'g.')
                        case 'DIR_All'
                            plot(X+ii,Y,'b.')
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
                        case 'DIR_All'
                            errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                            
                    end
                end
                
                
                % PLOT running mean ((all and stim catch) vs. (stim not catch)
                % plot no-stim
                X=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.timevals_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.([statfield '_SmMean']);
                
                plot(X+ii,Y,'k-','LineWidth',2);
                
                % plot stim
                if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
                    
                    X=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.timevals_sm;
                    X=X./24; % from hours to fraction of day.
                    Y=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.([statfield '_SmMean']);
                    
                    plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
                end
                
                % plot DIR
                if sum(strcmp(trialtypes,'DIR_All'))>0; % then this day has dir.
                    
                    X=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.timevals_sm;
                    X=X./24; % from hours to fraction of day.
                    Y=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.([statfield '_SmMean']);
                    
                    plot(X+ii,Y,'-','Color','b','LineWidth',2);
                end
                
                
                
                
                % plot mean of the day (if there are multipel trial types)
                if numtrialtypes>1;
                    Ytot=[];
                    for iii=1:numtrialtypes;
                        trialfield=trialtypes{iii};
                        Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                        Ytot=[Ytot Y];
                    end
                    
                    ffmean=mean(Ytot);
                    ffstd=std(Ytot);
                    ffn=length(Ytot);
                    ffsem=ffstd/sqrt(ffn-1);
                    
                    errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
                    
                end
                
                % Note down date
                Ylim=ylim;
                text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
                
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
            
            if ~isempty(DATSTRUCT.data{ii})
                
                % how many trial types today? e.g. StimCatch
                trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
                numtrialtypes=length(trialtypes);
                
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                    
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
                        case 'DIR_All'
                            errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                            
                    end
                end
                
                % plot mean of the day (if there are multipel trial types)
                if numtrialtypes>1;
                    Ytot=[];
                    for iii=1:numtrialtypes;
                        trialfield=trialtypes{iii};
                        Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                        Ytot=[Ytot Y];
                    end
                    
                    ffmean=mean(Ytot);
                    ffstd=std(Ytot);
                    ffn=length(Ytot);
                    ffsem=ffstd/sqrt(ffn-1);
                    
                    errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
                    
                end
                
                % Note down date
                Ylim=ylim;
                text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
                
            end
        end
    end
end

%% PLOT ALL DATA (SUBTRACTING BASELINE HOUR BY HOUR FLUCTUATION)
% DO THE SAME, BUT FOR BASELINE SUBTRACTED SONG --------------------
% PLOT - all vals, along with running mean (combining stimcatch and
% baseline) (at median time points).


% ==== PLOT WITH STIM CATCH AND NOT STIM COMBINED INTO ONE RUNNING AVERAGE
if (0)
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on
    title([timefield ' - ' statfield ' - (baseline subtracted, all renditions)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                
                
                switch trialfield
                    case 'All'
                        plot(X+ii,Y,'k.')
                    case 'StimCatch'
                        plot(X+ii,Y,'r.')
                    case 'StimNotCatch'
                        plot(X+ii,Y,'g.')
                    case 'DIR_All'
                        plot(X+ii,Y,'b.')
                        
                end
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                    case 'DIR_All'
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                end
            end
            
            
            % PLOT RUNGING mean ((all and stim catch) vs. (stim not catch)
            % plot no-stim
            X=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST_sm;
            X=X./24; % from hours to fraction of day.
            Y=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
            Ystd=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
            Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
            
            shadedErrorBar(X+ii,Y,Ysem,{'k-','LineWidth',2},1);
            
            %             plot(X+ii,Y,'k-','LineWidth',2);
            
            % plot stim
            if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                shadedErrorBar(X+ii,Y,Ysem,{'-','Color',[0.1 0.6 0.1],'LineWidth',2},1);
                
                %                 plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
            end
            
            % plot dir
            if sum(strcmp(trialtypes,'DIR_All'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                shadedErrorBar(X+ii,Y,Ysem,{'b-','LineWidth',2},1);
                
                %             plot(X+ii,Y,'-','Color','b','LineWidth',2);
            end
            
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    end
end


% == PLOT without combining stimcatch + baseline.
% baseline) (at median time points).
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title([timefield ' - ' statfield ' - (baseline subtracted, all renditions)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                
                
                switch trialfield
                    case 'All'
                        plot(X+ii,Y,'k.')
                    case 'StimCatch'
                        plot(X+ii,Y,'r.')
                    case 'StimNotCatch'
                        plot(X+ii,Y,'.', 'Color', 'b')
                    case 'DIR_All'
                        plot(X+ii,Y,'m.')
                        
                end
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
            end
            
            
            % PLOT RUNGING mean
            % plot stim off
            if isfield(DATSTRUCT.data{ii}.timewindow{i}, 'All');
                
                
                X=DATSTRUCT.data{ii}.timewindow{i}.All.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.All.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.All.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                shadedErrorBar(X+ii,Y,Ysem,{'k-','LineWidth',2},1);
            end
            
            % plot stim catch
            if any(strcmp(trialtypes,'StimCatch')); % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'r-','LineWidth',2},1);
                end                %             plot(X+ii,Y,'-','Color','b','LineWidth',2);
            end
            
            
            % plot stim
            if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'-','Color','b','LineWidth',2},1);
                end
                %                 plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
            end
            
            % plot dir
            if sum(strcmp(trialtypes,'DIR_All'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.timevals_fudgeDST_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.([statfield '_fudgeDST_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                if length(X)>1
                shadedErrorBar(X+ii,Y,Ysem,{'m-','LineWidth',2},1);
                end
                %             plot(X+ii,Y,'-','Color','b','LineWidth',2);
            end
            
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    lt_plot_annotation(1, 'red = catch', 'r');
end


% ==  PLOT ONLY MEANS
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title([timefield ' - Mean ' statfield ' (Baseline subtracted)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
            end
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    lt_plot_annotation(1, 'red = catch', 'r');
end



%% PLOT ALL DATA (WITHOUT SUBTRACTING BASELINE HOUR BY HOUR FLUCTUATION

% == PLOT without combining stimcatch + baseline.
% baseline) (at median time points).
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title([timefield ' - ' statfield ' - (not baseline subtr)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
%                 X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST;
%                 X=X./24; % from hours to fraction of day.
%                 Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                
                X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                
                
                switch trialfield
                    case 'All'
                        plot(X+ii,Y,'k.')
                    case 'StimCatch'
                        plot(X+ii,Y,'r.')
                    case 'StimNotCatch'
                        plot(X+ii,Y,'b.')
                    case 'DIR_All'
                        plot(X+ii,Y,'m.')
                        
                end
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
            end
            
            
            % ============================ PLOT RUNGING mean
            % plot stim off
            if isfield(DATSTRUCT.data{ii}.timewindow{i}, 'All');
                
                
                X=DATSTRUCT.data{ii}.timewindow{i}.All.timevals_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.All.([statfield '_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.All.([statfield '_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'k-','LineWidth',2},1);
                end
            end
            
            % plot stim catch
            if any(strcmp(trialtypes,'StimCatch')); % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.timevals_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.([statfield '_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.([statfield '_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'r-','LineWidth',2},1);
                                end
                %             plot(X+ii,Y,'-','Color','b','LineWidth',2);
            end
            
            
            % plot stim
            if sum(strcmp(trialtypes,'StimNotCatch'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.timevals_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.([statfield '_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.([statfield '_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'-','Color','b','LineWidth',2},1);
                                end                
                %                 plot(X+ii,Y,'-','Color',[0.1 0.6 0.1],'LineWidth',2);
            end
            
            % plot dir
            if sum(strcmp(trialtypes,'DIR_All'))>0; % then this day has stim.
                
                X=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.timevals_sm;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.([statfield '_SmMean']);
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.DIR_All.([statfield '_SmSTD']);
                Ysem=Ystd./sqrt(PARAMS.global.RunBin-1);
                
                                if length(X)>1;
                shadedErrorBar(X+ii,Y,Ysem,{'m-','LineWidth',2},1);
                                end                
                %             plot(X+ii,Y,'-','Color','b','LineWidth',2);
            end
            
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    lt_plot_annotation(1, 'red = catch', 'r');
end


% ==  PLOT ONLY MEANS
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title([timefield ' - Mean ' statfield ' (not baseline subtr)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
            end
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    lt_plot_annotation(1, 'red = catch', 'r');
end


%% SUBTRACT BASELINE EFFECT FOR EACH CATEGORY OF DATA (I.E. STIMNOTCATCH, STIMCATCH, NOTSTIM)



for i=1:length(TimeFieldsOfInterest); % for each window of interest

baseline_vals.All=[];
baseline_vals.StimCatch=[];
baseline_vals.StimNotCatch=[];
baseline_days=PARAMS.global.BaselineDays;

    for ii=baseline_days;
        
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
    
                rawvals=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield);
                
                % -- put into output                
                baseline_vals.(trialfield)=[baseline_vals.(trialfield) rawvals];
                
            end
        end
    end
    
    
    % ======= all baseline vals collected. get mean and subtract this mean from
    % all days
    
    
    % 1) mean of baselines
baseline_mean.All=mean(baseline_vals.All);
baseline_mean.StimCatch=mean(baseline_vals.StimCatch);
baseline_mean.StimNotCatch=mean(baseline_vals.StimNotCatch);
    
DATSTRUCT.baseline_data.timewindow{i}.All.([statfield])=baseline_vals.All;
DATSTRUCT.baseline_data.timewindow{i}.StimCatch.([statfield])=baseline_vals.StimCatch;
DATSTRUCT.baseline_data.timewindow{i}.StimNotCatch.([statfield])=baseline_vals.StimNotCatch;

    for ii=1:NumDays;
    
        if ~isempty(DATSTRUCT.data{ii})
           
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
 
               DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_MinBase'])=...
                   DATSTRUCT.data{ii}.timewindow{i}.(trialfield).(statfield)-baseline_mean.(trialfield);
            end
        end
    end
end


%% ====== PLOT ALL DATA (SUBTRACTING BASELINE SPECIFIC TO EACH TRIAL TYPE);
% note: skips any timewindow/trialtype conbination that lacks baseline data
% (just won't plot those data)

% == PLOT without combining stimcatch + baseline.
% baseline) (at median time points).
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    
    lt_figure; hold on; grid on;
    title([timefield ' - ' statfield ' - (each trial type subtr its own baseline)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
%                 X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.timevals_fudgeDST;
%                 X=X./24; % from hours to fraction of day.
%                 Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST']);
                


    if isempty(DATSTRUCT.baseline_data.timewindow{i}.All.([statfield]));
        continue;
    end


                X=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).timevals;
                X=X./24; % from hours to fraction of day.
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_MinBase']);
                
                
                switch trialfield
                    case 'All'
                        plot(X+ii,Y,'k.')
                    case 'StimCatch'
                        plot(X+ii,Y,'r.')
                    case 'StimNotCatch'
                        plot(X+ii,Y,'b.')
                    case 'DIR_All'
                        plot(X+ii,Y,'m.')
                        
                end
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+1+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
                
                
                % ==== PLOT RUNNING MEAN
                tmp=lt_running_stats(X, PARAMS.global.RunBin);
                Xsm=tmp.Median;
                
                tmp=lt_running_stats(Y, PARAMS.global.RunBin);
                Ysm=tmp.Mean;
                Ysm_sem=tmp.SEM;
                
                switch trialfield;
                    case 'All'
                        color='k';
                    case 'StimCatch'
                        color='r';
                    case 'StimNotCatch'
                        color='b';
                    case 'DIR_All'
                        color='m';
                end
                
                if length(Xsm)>1 & ~isnan(Ysm);
                shadedErrorBar(Xsm+ii,Ysm,Ysm_sem,{[color '-'],'LineWidth',2},1);
                end
                
            end
            
            
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_MinBase']);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+1+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
        lt_plot_annotation(1, 'red = catch', 'r');

end


% ==  PLOT ONLY MEANS
for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title([timefield ' - Mean ' statfield ' (each trial type subtr its own baseline)']);
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_MinBase']);
                
                % Plot mean and sem of this data type
                Y=Y(~isnan(Y));
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
                        errorbar(ii+iii*0.05,ffmean,ffsem,'bo','MarkerSize',9,'MarkerFaceColor','b')
                    case 'DIR_All'
                        errorbar(ii+iii*0.05,ffmean,ffsem,'mo','MarkerSize',9,'MarkerFaceColor','m')
                end
            end
            
            % plot mean of the day (if there are multipel trial types)
            if numtrialtypes>1;
                Ytot=[];
                for iii=1:numtrialtypes;
                    trialfield=trialtypes{iii};
                    Y=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).([statfield '_MinBase']);
                    Ytot=[Ytot Y];
                end
                
                Ytot=Ytot(~isnan(Ytot));
                ffmean=mean(Ytot);
                ffstd=std(Ytot);
                ffn=length(Ytot);
                ffsem=ffstd/sqrt(ffn-1);
                
                errorbar(ii+iii*0.05,ffmean,ffsem,'ko','MarkerSize',12)
                
            end
            
            % Note down date
            Ylim=ylim;
            text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
            
        end
    end
    lt_plot_annotation(1, 'red = catch', 'r');
end




%% PLOT FOR STIM ASSOC EXPERIMENTS (each day plot 2 datapoints of means)
if exist('MetaParams','var');
    if isfield(MetaParams,'AssocExpt');
        
        ExperimentTypes=fieldnames(MetaParams.AssocExpt.Timeline); % e.g. Stim + WNup
        
        for i=1:length(TimeFieldsOfInterest); % for each window of interest
            timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
            
            figure; hold on;
            
            % plot each day
            Ymin=[]; % collect min and max
            Ymax=[];
            
            for ii=1:NumDays;
                
                % skip baseline days
                if ~any(BaselineDays==ii)
                    
                    % what's today's date?
                    todaydatenum=PARAMS.global.FirstDate_num+ii-1;
                    todaydate_datestr=datestr(todaydatenum,'ddmmmyyyy');
                    
                    
                    
                    % what kind of experiment was today?
                    experimenttype=[];
                    for kk=1:length(ExperimentTypes);
                        if any(strcmp(MetaParams.AssocExpt.Timeline.(ExperimentTypes{kk}),todaydate_datestr));
                            % found what experiment today was. save info and break
                            % loop.
                            experimenttype=ExperimentTypes{kk};
                            break
                        end
                    end
                    
                    % if today has no meta description, skip.
                    if isempty(experimenttype);
                        continue
                    end
                    
                    
                    % based on today;'s experiment, decide which subplot to
                    % use.
                    switch experimenttype
                        case {'Stim_WNup_UniDir', 'Stim_WNup_BiDir'}
                            subplotnumber=1;
                        case 'Stim_WNdn_UniDir'
                            subplotnumber=2;
                        case 'noStim_WNup'
                            subplotnumber=3;
                        case 'noStim_WNdn'
                            subplotnumber=4;
                        case 'Stim_WNup_probe'
                            subplotnumber=5;
                        case 'Stim_WNdn_probe'
                            subplotnumber=6;
                    end
                    hsplot(subplotnumber)=subplot(1,6,subplotnumber); hold on;
                    title(experimenttype);
                    
                    % PLOT
                    Ystim=[];
                    Ynostim=[];
                    Yall=[];
                    if ~isempty(DATSTRUCT.data{ii})
                        
                        Ystim=DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.MINUSBaseHrBins.([statfield '_fudgeDST']);
                        Ynostim=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.MINUSBaseHrBins.([statfield '_fudgeDST']);
                        
                        % Plot mean and sem of this data type
                        Ystim=Ystim(~isnan(Ystim));
                        Ystim_mean=mean(Ystim);
                        
                        Ynostim=Ynostim(~isnan(Ynostim));
                        Ynostim_mean=mean(Ynostim);
                        
                        Yall=[Ynostim_mean Ystim_mean];
                        
                        plot([1 2], Yall,'-k');
                        
                        plot(2,Yall(2),'go','MarkerSize',9,'MarkerFaceColor','g');
                        plot(1,Yall(1),'ro','MarkerSize',9,'MarkerFaceColor','r');
                        
                        Ymax=max([Ymax Yall]);
                        Ymin=min([Ymin Yall]);
                        
                        
                        %             % Note down date
                        %             Ylim=ylim;
                        %             text(ii+0.5, Ylim(2), datestr(FirstDate_num+ii-1,'ddmmmyyyy'))
                        
                    end
                end
            end
            
            % set all subplots to have same max and min axes
            for jjj=1:6;
                subplot(1,6,jjj);
                ylim([Ymin Ymax]);
            end
            %             linkaxes(hsplot,'y')
            subtitle([timefield ' - Mean ' statfield ' (Baseline subtracted)']);
        end
        
    end
end



%% PLOT RUNNING VARIABILITY MEAN (PER DAY)

for i=1:length(TimeFieldsOfInterest); % for each window of interest
    timefield=PARAMS.global.TimeFields.char{TimeFieldsOfInterest(i)};
    
    lt_figure; hold on; grid on;
    title(['mean of STD (over time bins); ' timefield ' - ' statfield ' - (baseline subtracted, all renditions)']);
    ylabel('mean std (hz)');
    xlabel('day');
    
    % plot each day
    for ii=1:NumDays;
        if ~isempty(DATSTRUCT.data{ii})
            
            
            % how many trial types today? e.g. StimCatch
            trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
            numtrialtypes=length(trialtypes);
            
            for iii=1:numtrialtypes;
                trialfield=trialtypes{iii};
                
                Ystd=DATSTRUCT.data{ii}.timewindow{i}.(trialfield).MINUSBaseHrBins.([statfield '_fudgeDST_SmSTD']);
                Ystd=Ystd(~isnan(Ystd));
                meanYstd=mean(Ystd);
                
                switch trialfield
                    case 'All'
                        plot(ii,meanYstd,'ko','MarkerSize',9,'MarkerFaceColor','k')
                    case 'StimCatch'
                        plot(ii,meanYstd,'ro','MarkerSize',9,'MarkerFaceColor','r')
                    case 'StimNotCatch'
                        plot(ii,meanYstd,'go','MarkerSize',9,'MarkerFaceColor','g')
                    case 'DIR_All'
                        plot(ii,meanYstd,'bo','MarkerSize',9,'MarkerFaceColor','b')
                end
            end
        end
    end
end


%% TROUBLESHOOTING - check to make sure stim epochs extracted are correct

% figure(16);
% hold on;
%
% for ii=1:NumDays;
%     if ~isempty(StimEpochs{ii});
%
%         NumEpochs=length(StimEpochs{ii}.epoch);
%         for k=1:NumEpochs;
%
%             % plot stim
%             SCinds=StimEpochs{ii}.epoch{k}.StimCatch;
%             Tvals=DATSTRUCT.data{ii}.timewindow{i}.StimCatch.timevals;
%
%             plot(ii+Tvals(SCinds)./24,200+10*rand,'xr');
%
%             % plot stim pre
%             try
%             SCinds=StimEpochs{ii}.epoch{k}.All_preceding;
%             Tvals=DATSTRUCT.data{ii}.timewindow{i}.All.timevals;
%
%             plot(ii+Tvals(SCinds)./24,200+10*rand,'xk');
%             catch err
%             end
%
%             % plot stim post
%             try
%             SCinds=StimEpochs{ii}.epoch{k}.All_post;
%             Tvals=DATSTRUCT.data{ii}.timewindow{i}.All.timevals;
%
%             plot(ii+Tvals(SCinds)./24,200+10*rand,'xk');
%             catch err
%
%             end
%
%
%
%             disp(['day: ' num2str(ii)]);
%             disp(['epoch: ' num2str(k)]);
%
%             pause
%
%         end
%     end
% end


%% FOR ALL STIM EPOCHS, plot reversion versus various things

if plotStimEpochs==1;
    
    
    % == PLOT STIM EPOCHS, ETC.
    RunBin=4;
    NumEdgeTrials=4;
    StimEpochs_aligned=lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs1(PARAMS, DATSTRUCT, StimEpochs, RunBin, NumEdgeTrials, TimeFieldsOfInterest, statfield);
    
    % == REGRESSIONS, day level stats
    NumStimRends=25; % i.e. 25 stim and 25 nonstim.
    [StimEpochs_aligned, Y_day_level_regression]=lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs2(StimEpochs_aligned, PARAMS, DATSTRUCT, NumStimRends, TimeFieldsOfInterest, statfield);
    
    % == REGRESSIONS, one datapoint for each stim epoch (and stats based on
    % that epoch)
    % NOTE: GO IN HERE TO PERFORM PERMUTATION TESTS
    NumTrials_prestim=200;
    NumTrials_stim=50;
    Conditions.Only_Plot_Days_With_One_Stim=0;
    Conditions.OnlyPlotFirstEpochOfDay=0;
    Conditions.OnlyPlotIfNoStimYesterday=0;
    
    [StimEpochs_aligned_FINAL, RealStats] = lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs3(StimEpochs_aligned, StimEpochs, DATSTRUCT, NumTrials_prestim, NumTrials_stim, TimeFieldsOfInterest, statfield, Conditions);
    
    
    
    % == OLD, IGNORE
    % BinSize=15; % for running avg
    % NumTrials_reversion=10; % e.g. to calculate reversion of stim catch get 20. (make empty if want all)
    % NumTrials_otherstats=50; % trials to take for things regressing against reversion (e.g. slope of learning etc).
    % NumTrials_plot=150; % trials to plot in scatter plot
    
    % lt_Opto_Stim_analy_SUMMARY_PlotOverTime_sub1;
    % === STOP IGNORING
    
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % ========== PERFORM REGRESSIONS FOR VARYING SIZE WINDOWS
    Conditions.Only_Plot_Days_With_One_Stim=1;
    Conditions.OnlyPlotFirstEpochOfDay=1;
    Conditions.OnlyPlotIfNoStimYesterday=1;
    
    PreStim_List=[5 10 20 40 60 80 100 200 500];
    DurStim_List=[5 10 20 40 60 100];
    % PreStim_List=[10 20 30 40];
    % DurStim_List=[10 20];
    
    % make matrix of all the combinations of Pre and dur stim
    clear Collected_regression_stats;
    
    for i=1:length(PreStim_List);
        for ii=1:length(DurStim_List);
            
            % close figures
            close all
            
            % === perform regression using corresponding Trial numbers
            NumTrials_prestim=PreStim_List(i);
            NumTrials_stim=DurStim_List(ii);
            
            [~, RealStats] = lt_Opto_Stim_analy_SUMMARY_PlotOverTime_StimEpochs3(StimEpochs_aligned, StimEpochs, DATSTRUCT, NumTrials_prestim, NumTrials_stim, TimeFieldsOfInterest, statfield, Conditions);
            
            % === collect stats for each time window
            for iii=1:length(TimeFieldsOfInterest);
                Collected_regression_stats{iii}.stats(i,ii) = RealStats.timewindow{iii}.regressions;
            end
            
        end
    end
    
    
    % ====== PLOT COLLECTED REGRESSION STATS
    RegressionStat_list=fieldnames(Collected_regression_stats{1}.stats(1,1));
    Pvals=[];
    R2vals=[];
    Slopevals=[];
    P_less_than_alpha=[];
    
    for jj=1:length(RegressionStat_list);
        regressionstat=RegressionStat_list{jj};
        %     regressionstat='recent_learning_StimCombined';
        
        % === PLOT P VALUE
        lt_figure; hold on;
        for iii=1:length(TimeFieldsOfInterest);
            lt_subplot(ceil(length(TimeFieldsOfInterest)/2), 2, iii); hold on;
            title(['TimeWindow: ' num2str(TimeFieldsOfInterest(iii))]);
            
            % == Collect p and R2 values across all bin sizes
            for i=1:length(PreStim_List);
                for ii=1:length(DurStim_List);
                    
                    Pvals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).p;
                    R2vals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).r2;
                    Slopevals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).slope;
                    
                    % -- find any with p<0.05
                    if Collected_regression_stats{iii}.stats(i,ii).(regressionstat).p<0.05;
                        P_less_than_alpha(i,ii)=1;
                    else
                        P_less_than_alpha(i,ii)=0;
                    end
                    
                end
            end
            
            % == plot 3d
            imagesc(Pvals,[0 0.15]);
            % %    plot(P_less_than_alpha, 'x', 'MarkerSize', 10)
            %       h=imagesc(P_less_than_alpha);
            colormap('jet')
            colorbar
            
            % == annotate
            set(gca, 'XTick', 1:length(DurStim_List));
            set(gca,'XTickLabel', DurStim_List);
            
            set(gca, 'YTick', 1:length(PreStim_List));
            set(gca,'YTickLabel', PreStim_List);
            %
            %             xlim([DurStim_List(1) DurStim_List(end)]);
            %             ylim([PreStim_List(1) PreStim_List(end)]);
            %
            
            ylabel('Pre-Stim bin size');
            xlabel('Dur stim bin size');
        end
        lt_subtitle(['p-value: ' regressionstat])
        
        % ===== PLOT SLOPE
        lt_figure; hold on;
        for iii=1:length(TimeFieldsOfInterest);
            lt_subplot(ceil(length(TimeFieldsOfInterest)/2), 2, iii); hold on;
            title(['TimeWindow: ' num2str(TimeFieldsOfInterest(iii))]);
            
            % == Collect p and R2 values across all bin sizes
            for i=1:length(PreStim_List);
                for ii=1:length(DurStim_List);
                    
                    Pvals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).p;
                    R2vals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).r2;
                    Slopevals(i,ii)=Collected_regression_stats{iii}.stats(i,ii).(regressionstat).slope;
                    
                    % -- find any with p<0.05
                    if Collected_regression_stats{iii}.stats(i,ii).(regressionstat).p<0.05;
                        P_less_than_alpha(i,ii)=1;
                    else
                        P_less_than_alpha(i,ii)=0;
                    end
                end
            end
            
            % == plot 3d
            imagesc(Slopevals);
            % %    plot(P_less_than_alpha, 'x', 'MarkerSize', 10)
            %       h=imagesc(P_less_than_alpha);
            colormap('jet')
            colorbar
            
            
            % == annotate
            set(gca, 'XTick', 1:length(DurStim_List));
            set(gca,'XTickLabel', DurStim_List);
            
            set(gca, 'YTick', 1:length(PreStim_List));
            set(gca,'YTickLabel', PreStim_List);
            
            %             xlim([DurStim_List(1) DurStim_List(end)]);
            %             ylim([PreStim_List(1) PreStim_List(end)]);
            
            
            ylabel('Pre-Stim bin size');
            xlabel('Dur stim bin size');
        end
        lt_subtitle(['Slope: ' regressionstat])
        
    end
end


%% TO DO : Plot stim aligned, trial by trial (i.e. contours)

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



%% below was for baseline subtraacted IGNORE
if (0)
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
    
end

%% DOES REVERSION CORRELATE WITH 1) amount of learning 2) slope preceding, 3) slope during stim epoch

if (0)
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
            
            if ~isempty(DATSTRUCT.data{ii})
                
                if ~any(PARAMS.global.BaselineDays==ii); % ignore baseline days
                    % how many trial types today? e.g. StimCatch
                    trialtypes=fieldnames(DATSTRUCT.data{ii}.timewindow{i});
                    
                    % only continue if today has both stim catch and notcatch
                    if any(strcmp(trialtypes,'StimCatch')) && any(strcmp(trialtypes,'StimNotCatch'));
                        
                        % calculate reversion
                        Ystim=mean(DATSTRUCT.data{ii}.timewindow{i}.StimNotCatch.([statfield]));
                        Ycatch=mean(DATSTRUCT.data{ii}.timewindow{i}.StimCatch.([statfield]));
                        
                        REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.([statfield]) = Ystim-Ycatch;
                        
                        % collect all data
                        ReversionTot=[ReversionTot REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.([statfield])];
                        
                        % PLOT
                        % 1) Reversion vs. learning in that day
                        learning = mean(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])); % learning as mean of vals subtracted baseline
                        
                        figure(hfig1)
                        plot(learning,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield),'ko','MarkerFaceColor','k');
                        text(learning+1,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                        
                        % 2) change in learning (today minus yesterday);
                        if ii~=1 && ~isempty(DATSTRUCT.data{ii-1})
                            
                            learningDiff=mean(DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']))...
                                -mean(DATSTRUCT_2.data{ii-1}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST']));
                            LearningDiffTot=[LearningDiffTot learningDiff];
                            ReversionForLearningDiff=[ReversionForLearningDiff REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield)];
                            
                            figure(hfig2);
                            plot(learningDiff,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield),'ro','MarkerFaceColor','r');
                            text(learningDiff+1,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                        end
                        
                        % 3) Within day slope (all non-stim trials).
                        Tvals=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.timevals_fudgeDST;
                        x=[ones(length(Tvals),1) Tvals']; % versus time
                        y=DATSTRUCT_2.data{ii}.timewindow{i}.All_StimCatch_combined.MINUSBaseHrBins.([statfield '_fudgeDST'])';
                        [b,~,~,~,~]=regress(y,x); % get slope (b2)
                        
                        slope=b(2);
                        
                        figure(hfig3);
                        plot(slope,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield),'go','MarkerFaceColor','g');
                        text(slope+1,REVERSION.data{ii}.timewindow{i}.DayMeans.StimMinusCatch.(statfield)+1,[num2str(ii)])
                        
                        
                        % collect all data
                        LearningTot=[LearningTot learning];
                        WithinDaySlopeTot=[WithinDaySlopeTot slope];
                        
                        
                    end
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
    
end

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





%% SAVE

PARAMS.global.savedir=[PARAMS.global.BirdDir '/Opto_Stim_analy'];

try cd(PARAMS.global.savedir);
catch err
    mkdir(PARAMS.global.savedir);
    cd(PARAMS.global.savedir);
end

% save in subdir based on experiment.
uscores=strfind(PARAMS.global.ListOfDirs{1}, '_');
PARAMS.global.phrase=PARAMS.global.ListOfDirs{1}(uscores(1)+1:uscores(2)-1); % phrase, based on dir name of first dir
try cd(PARAMS.global.phrase);
catch err
    mkdir(PARAMS.global.phrase);
    cd(PARAMS.global.phrase);
end

% Save structures
save('PARAMS','PARAMS');
save('DATSTRUCT','DATSTRUCT');
save('DATSTRUCT_2','DATSTRUCT_2');
save('StimEpochs', 'StimEpochs');

% Make note
tstamp=lt_get_timestamp(0);
fid=fopen(['DONE_' tstamp '.txt'],'w');
fclose(fid);


if plotStimEpochs==1;
    StimEpochs_aligned=StimEpochs_aligned_FINAL;
    save('StimEpochs_aligned','StimEpochs_aligned');
    save('Y_day_level_regression','Y_day_level_regression');
end

% save FIgures
try cd('FIGURES');
catch err
    mkdir('FIGURES');
    cd('FIGURES');
end

tstamp=lt_get_timestamp(0);
try cd(tstamp);
catch err
    mkdir(tstamp);
    cd(tstamp);
end

lt_save_all_figs;

cd('../../')
disp('DONE AND SAVED!');

    
    





