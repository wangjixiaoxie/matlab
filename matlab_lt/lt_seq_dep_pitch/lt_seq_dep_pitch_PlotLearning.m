function [Params, AllDays_PlotLearning, AllDays_RawDatStruct]=lt_seq_dep_pitch_PlotLearning(Params, AllDays_RawDatStruct,saveON)
%% 5/26/15 - now first checks if data can be separated into muscimol vs. PBS.  if latter is true, then first separates data before running normal plot learning stuff
% does this for both musc and PBS data.

% checks by looking for entries in Params.PlotLearning.MuscimolSchedule

%% LT 5/19/15 - Added stuff to analyze muscimol epochs. 
% UPDATED PARAMS entry (showing all):

% Params.PlotLearning.MuscimolSchedule={...
%     {'03May2015', '1337', '1745'}, ...
%     {'04May2015', '1259', '1700'}, ...
%     {'05May2015', '1401', '2036'}, ...
%     {'07May2015', '1125', '1719'}, ...
%     {'08May2015', '1155', '1636'}, ...
%     {'09May2015', '1404', '1915'}, ...
%     {'10May2015', '1238', '1831'}, ...
%     {'11May2015', '1215', '1555'}, ...
%     {'18May2015', '1150', '1643'}};

% Params.PlotLearning.MuscimolDaysToThrowout={'08May2015', '10May2015'};

% Params.PlotLearning.plotWNdays=1; % if 1, then plots WN lines, if 0, then no.
% Params.PlotLearning.DayBinSize=3; % 3 day running avg.
% saveON=1;

% NOTE: saves muscimol information into raw dat and plot learning structs,
% and saves both of those to overwrite old




%% LT 4/28/15 - copied to "OLD" and here continuing editing - 
% now saves in all one directiory, and overwrites old stuff, and no
% timstamp in name.


%% LT 2/6/15 - Given raw dat struct compiled over days (from lt_seq_dep_pitch_DayRawDat) plots learning, and outputs stats struct
% Converted from lt_compile_seq_dep_pitch_data_PLOTDirUndir.
% DIFFERENCE: here plots one structure (e.g. undirected song). Old code
% also compared to directed song, so would load two structures.  Can easily
% plot to compare after this code outputs stats structure



%% PRE-PROCESSING

% ==== DEFAULTS
if ~exist('saveON','var');
    saveON=1;
end

save_AllDays_RawDat=0; % will save only if this is the first time running this code and if performed parsing of PBS/MUSC.

if ~isfield(Params.SeqFilter.SylLists, 'FieldsToPlot');
    Params.SeqFilter.SylLists.FieldsToPlot{1}=[Params.SeqFilter.SylLists.SylsSame Params.SeqFilter.SylLists.SylsDifferent];
end


% Extract params
NumDays=length(AllDays_RawDatStruct); % total days
SylFieldsAll=fieldnames(AllDays_RawDatStruct{1}.data); % all seq and syls (assume 1st day has all possible syls)
FirstDay=Params.SeqFilter.FirstDay;
LastDay=Params.SeqFilter.LastDay;
plotWNdays=Params.PlotLearning.plotWNdays;




% Check other days to make sure SylFieldsAll captures all syls for all days
for i=1:NumDays;
    if ~isempty(AllDays_RawDatStruct{i})
        tmp=fieldnames(AllDays_RawDatStruct{i}.data);
        if length(tmp)>length(SylFieldsAll); % i.e. this day has more syls
            disp(['PROBLEM - day ' num2str(i) ' has more syls than baseline day 1 - replacing with new SylFieldsAll (not a problem)']);
            % replace SylFieldsAll with updated version with more syls
            SylFieldsAll=tmp;
        elseif length(tmp)<length(SylFieldsAll);
            disp(['NOTE - day ' num2str(i) ' has fewer syls than baseline day 1 (not a problem)']);
            disp('that day:');
            disp(tmp)
            disp('baseline:');
            disp(SylFieldsAll);
        end
    end
end


% For each syl, is is single syl or sequence?
X=[];
for i=1:length(SylFieldsAll);
    X(i)=length(SylFieldsAll{i});
end
SylFieldsSingle=SylFieldsAll(X==1); % only single syls
SylFieldsSeq=SylFieldsAll(X>1); % only sequences


% =============== Get relevant syls - i.e. only syls uniquely defined by sequence context
% (i.e. don't double count b and a[b]);
SylFields_Unique={};
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder)
    SylFields_Unique=[SylFields_Unique Params.SeqFilter.SylLists.FieldsInOrder{i}];
end

% some syls are unique but not specified in sequence.  extract those as
% well.

UpperSyls=regexp(SylFields_Unique, '[A-Z]'); % find all upper case syls alrady in unique syls. if any single syls are not represented, then put the single syl into unique syls
upper_syls_lowered={};
for i=1:length(UpperSyls);
    if isempty(UpperSyls{i});
        upper_syls_lowered{i}=SylFields_Unique{i};
    else
        upper_syls_lowered{i}=lower(SylFields_Unique{i}(UpperSyls{i}));
    end
end

% for all single syls, if not represented in lower syls (i.e. syls that
% are already called unqiue), then add it to the unique group
for i=1:length(SylFieldsSingle);
    % only continue if the single syl is lowercase (upper is redundant
    % as has to have a lower case of same thing)
    if regexp(SylFieldsSingle{i}, '[A-K]');
        continue
    end
    
    if ~any(strcmp(SylFieldsSingle{i}, upper_syls_lowered));
        SylFields_Unique=[SylFields_Unique SylFieldsSingle{i}];
    end
end
% ====================================================




% =========== Convert WNdays to index days
global WNTimeOnInd % make global, so can use in subfunction below.
global WNTimeOffInd

if Params.PlotLearning.plotWNdays==1;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeON});
    WNTimeOnInd=X.JustDays_rel;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeOFF});
    WNTimeOffInd=X.JustDays_rel;
end


% ===== Get baseline days
if isfield(Params.SeqFilter, 'BaselineDays');
    disp(['replacing Params.SeqFilter.BaselineDays from ' num2str(Params.SeqFilter.BaselineDays) ' to ' ...
        num2str(1:WNTimeOnInd-1)]);
end
Params.SeqFilter.BaselineDays=1:WNTimeOnInd-1;



% ========= Get days to mark, if exist.
global DaysToMarkInds
DaysToMarkInds={};
X={};
isfield(Params.SeqFilter,'DaysToMark')
if isfield(Params.SeqFilter,'DaysToMark'); % if there are specific days to average over and look at.
    DaysToMark=Params.SeqFilter.DaysToMark;
    for i=1:length(DaysToMark);
        X{i}=lt_convert_EventTimes_to_RelTimes(FirstDay,{DaysToMark{i}});
        DaysToMarkInds{i}=X{i}.JustDays_rel;
    end
end


% For saving
timestampSv=lt_get_timestamp(0);
SaveDir=Params.SeqFilter.savedir;


% Are there muscimol days?
plotLMANmusc=0;
if isfield(Params.PlotLearning, 'MuscimolSchedule');
    % then there might be muscimol days
    
    if ~isempty(Params.PlotLearning.MuscimolSchedule);
        plotLMANmusc=1;
    end
end

    
% IF THERE ARE MUSCIMOL DAYS, FIGURE OUT STARTS AND STOPS
if plotLMANmusc==1;
    
    % === WHAT ARE SWITCH TIMES FOR EACH DAY?
   for i=1:length(Params.PlotLearning.MuscimolSchedule);

       
    date_today=Params.PlotLearning.MuscimolSchedule{i}{1};
    start_time=Params.PlotLearning.MuscimolSchedule{i}{2};
    end_time=Params.PlotLearning.MuscimolSchedule{i}{3};
    
    % convert date to ind
    event_times=lt_convert_EventTimes_to_RelTimes(FirstDay, {date_today});
    date_ind=event_times.JustDays_rel;
    
    % throw out if less than 1
    if date_ind<1;
        continue
    end
    
    % convert times to hours
    [~, B] = lt_convert_datenum_to_hour(datenum(start_time,'HHMM'));
    Params.PlotLearning.MuscimolSchedule_ByDayInds{date_ind}.start=B.hours;
    
    [~, B] = lt_convert_datenum_to_hour(datenum(end_time,'HHMM'));
    Params.PlotLearning.MuscimolSchedule_ByDayInds{date_ind}.end=B.hours;
   end

   % make sure every day is either empty, or has ind
   if length(Params.PlotLearning.MuscimolSchedule_ByDayInds)<NumDays; % if not enough days, then likely because missing last day
       missing_days=length(Params.PlotLearning.MuscimolSchedule_ByDayInds)+1:NumDays;
       for l=missing_days;
           Params.PlotLearning.MuscimolSchedule_ByDayInds{l}=[];
       end
   end
   
   % == WHAT IS MEDIAN MUSCIMOL SWITCH TIME?
   Y=[];
   for i=1:length(Params.PlotLearning.MuscimolSchedule);
       try
           switch_time=Params.PlotLearning.MuscimolSchedule{i}{2};
           [~, B]=lt_convert_datenum_to_hour(datenum(    switch_time, 'HHMM'));
           
           Y=[Y B.hours];
           
           
       catch err % nothing today
       end
   end
   
   Params.PlotLearning.MuscimolSchedule_MedianStartTime=median(Y);
   
   % == WHat are days to throw out (ni.e. muscimol failed)
   Params.PlotLearning.MuscimolDaysToThrowout_Inds=[];
   if isfield(Params.PlotLearning, 'MuscimolDaysToThrowout');
       if ~isempty(Params.PlotLearning.MuscimolDaysToThrowout);
           tmp=lt_convert_EventTimes_to_RelTimes(FirstDay, Params.PlotLearning.MuscimolDaysToThrowout);
           Params.PlotLearning.MuscimolDaysToThrowout_Inds=tmp.JustDays_rel;
       end
   end
end


% UPDATE PARAMS STRUCTURE
Params.PlotLearning.SylFieldsAll=SylFieldsAll;
Params.PlotLearning.SylFieldsSingle=SylFieldsSingle;
Params.PlotLearning.SylFieldsSeq=SylFieldsSeq;
Params.PlotLearning.WNTimeOnInd=WNTimeOnInd;
Params.PlotLearning.WNTimeOffInd=WNTimeOffInd;
Params.PlotLearning.DaysToMarkInds=DaysToMarkInds;
Params.PlotLearning.savedir=SaveDir;
Params.PlotLearning.SylFields_Unique=SylFields_Unique;


%% ==== CHANGE FIELDS TO PLOT TO Unique Syls (similar vs different)
Params.SeqFilter.SylLists=rmfield(Params.SeqFilter.SylLists, 'FieldsToPlot');

Params.SeqFilter.SylLists.FieldsToPlot{1}=[Params.SeqFilter.SylLists.TargetSyls Params.SeqFilter.SylLists.SylsSame];
if isfield(Params.SeqFilter.SylLists, 'SylsDifferent');
    Params.SeqFilter.SylLists.FieldsToPlot{2}=[Params.SeqFilter.SylLists.SylsDifferent];
    combined_fields=[Params.SeqFilter.SylLists.FieldsToPlot{1} Params.SeqFilter.SylLists.FieldsToPlot{2}];
end
combined_fields=[Params.SeqFilter.SylLists.FieldsToPlot{1}];


% any additional things in Unique syls that is not there, put in syls
% different
tmp=setdiff(Params.PlotLearning.SylFields_Unique, combined_fields); % pulls out unique sysl that are not in one of fileds to plot
if ~isempty(tmp);
    if length(Params.SeqFilter.SylLists.FieldsToPlot)>1;
        Params.SeqFilter.SylLists.FieldsToPlot{2}=[Params.SeqFilter.SylLists.FieldsToPlot{2}, tmp];
    else
        Params.SeqFilter.SylLists.FieldsToPlot{2}=tmp;
    end
end

%% IF THERE ARE MUSCIMOL DATA, THEN SEPARATE PBS VS. MUSCIMOL DATA BEFORE ANALYZING
% PBS data will slot into the normal analyses.  will have additional field
% holding MUSC data that can also be analyzed using same plot learning
% code.

if plotLMANmusc==1;
    save_AllDays_RawDat=1;
    
    % ======================= 
    % Check if already done this analysis. If so, then do not do it again -
    % that leads to error.
    
    Parse_Musc_Pbs=1; % default
    if isfield(Params.PlotLearning, 'ExptCondition_codes');
        % then I know this has been done before;
        Parse_Musc_Pbs=0; % don't parse
        save_AllDays_RawDat=0; % don't bother resaving.
    end
    
    if Parse_Musc_Pbs==1;
        
        ExptCondition_codes={'PBS','MUSC'}; % i.e. if filename says "PBS", then is code 1 (musc=2 and so forth)
        Params.PlotLearning.ExptCondition_codes=ExptCondition_codes;
        
        for ii=1:NumDays;
            if isempty(AllDays_RawDatStruct{ii});
                % then no data today, skip.
                continue
            end
            
            
            % == 1) ANNOTATE each rendition with experimental condition.
            for j=1:length(SylFieldsAll);
                syl=SylFieldsAll{j};
                
                
                % ==== SORT TODAY's TRIALS INTO CONDITIONS ==================
                % get all filenames
                filenames=[];
                if isfield(AllDays_RawDatStruct{ii}.data, syl);
                    filenames=AllDays_RawDatStruct{ii}.data.(syl)(:,5);
                end
                
                if isempty(filenames);
                    continue;
                end
                
                % get conditions from filenames
                ExptConditions=[];
                for iii=1:length(filenames);
                    underscores=strfind(filenames{iii},'_');
                    exptcond=filenames{iii}(underscores(1)+1:underscores(2)-1);
                    
                    % convert that condition (string) to code (1, 2, ...);
                    if any(strcmp(ExptCondition_codes, exptcond));
                        ExptConditions(iii)=find(strcmp(ExptCondition_codes, exptcond));
                    else
                        disp('PROBLEM - a song is not PBS or MUSC based on filename');
                    end
                end
                
                
                
                % == SAVE CONDITION TO RAW DAT
                % STRUCTURES
                % 1) Raw Dat structure
                AllDays_RawDatStruct{ii}.data.(syl)(:,12)=num2cell(ExptConditions');
                
                % update legend in params
                Params.DayRawDat.Legend_data{1,:}={'1_FF','2_PitchContour','3_SoundDat','4_Spec',...
                    '5_FileName','6_DateNum','7_LabelStr','8_NotePos','9_Trig','10_NoteDur','11_CatchTrial','12_ExptCond'};
                Params.DayRawDat.ExptCondition_codes=ExptCondition_codes;
                
            end
            
            % == 2) backup all data for this dat into a new field, and erase the old
            % field
            AllDays_RawDatStruct{ii}.data_PBS_and_MUSC=AllDays_RawDatStruct{ii}.data; % backup
            AllDays_RawDatStruct{ii}=rmfield(AllDays_RawDatStruct{ii}, 'data'); % erase
            
            
            % == 3) Refill data structure, with just PBS data. (to separate
            % data based on PBS vs. MUSC)
            % ONLY KEEP UNIQUE SYLS
            for j=1:length(SylFields_Unique);
                syl=SylFields_Unique{j};
                
                % -- find inds for PBS data
                PBS_indicator=find(strcmp(ExptCondition_codes,'PBS'));
                
                if isfield(AllDays_RawDatStruct{ii}.data_PBS_and_MUSC, syl);
                    % then there are data
                    PBS_inds=cell2mat(AllDays_RawDatStruct{ii}.data_PBS_and_MUSC.(syl)(:,12))==PBS_indicator;
                    
                    % -- slot data of those inds back into "data" field
                    AllDays_RawDatStruct{ii}.data.(syl)=AllDays_RawDatStruct{ii}.data_PBS_and_MUSC.(syl)(PBS_inds,:);
                end
                % == Filter out Musc data
                % -- find inds for MUSC data
                MUSC_indicator=find(strcmp(ExptCondition_codes,'MUSC'));
                
                if isfield(AllDays_RawDatStruct{ii}.data_PBS_and_MUSC, syl)
                    MUSC_inds=cell2mat(AllDays_RawDatStruct{ii}.data_PBS_and_MUSC.(syl)(:,12))==MUSC_indicator;
                    
                    % -- slot data of those inds back into "data_MUSC" field
                    AllDays_RawDatStruct{ii}.data_MUSC.(syl)=AllDays_RawDatStruct{ii}.data_PBS_and_MUSC.(syl)(MUSC_inds,:);
                end
                
                
            end
            
            % == 4) REMOVE THE BACKUP - to save disk space during saving.
            AllDays_RawDatStruct{ii}=rmfield(AllDays_RawDatStruct{ii}, 'data_PBS_and_MUSC');
        end
    end
end




%% PROCESS DATA - extract means, etcs
if plotLMANmusc==1; % if LMAN data exists, then will do this for both PBS and MUSC
    [EpochData, DataMatrix, DataMatrix_Targ, Params] = lt_seq_dep_pitch_PlotLearning_PreProcess(Params, AllDays_RawDatStruct, SylFields_Unique); % PBS DATA
    
    [EpochData_MUSC, DataMatrix_MUSC, DataMatrix_Targ_MUSC, Params] = lt_seq_dep_pitch_PlotLearning_PreProcess(Params, AllDays_RawDatStruct, SylFields_Unique, 1); % MUSC DATA

% == PUT THOSE into AllDaysPlotLearning struct
AllDays_PlotLearning.EpochData=EpochData;
AllDays_PlotLearning.EpochData_MUSC=EpochData_MUSC;
AllDays_PlotLearning.DataMatrix=DataMatrix;
AllDays_PlotLearning.DataMatrix_MUSC=DataMatrix_MUSC;
AllDays_PlotLearning.DataMatrix_Targ=DataMatrix_Targ;
AllDays_PlotLearning.DataMatrix_Targ_MUSC=DataMatrix_Targ_MUSC;


else
    [EpochData, DataMatrix, DataMatrix_Targ, Params] = lt_seq_dep_pitch_PlotLearning_PreProcess(Params, AllDays_RawDatStruct, SylFields_Unique);

% == PUT THOSE into AllDaysPlotLearning struct
AllDays_PlotLearning.EpochData=EpochData;
AllDays_PlotLearning.DataMatrix=DataMatrix;
AllDays_PlotLearning.DataMatrix_Targ=DataMatrix_Targ;

end





%% PLOT - Day means over learning


if plotLMANmusc==1; % IF LMAN DATA EXISTS, THEN PLOT SEPARATE FOR MUSC AND PBS
    datamatrix_fields={'DataMatrix', 'DataMatrix_MUSC'};
else
    datamatrix_fields={'DataMatrix'};
end


for dd=1:length(datamatrix_fields);
    datamatfield=datamatrix_fields{dd};
    
    lt_figure; hold on;
    
    for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        % Plot figure for this set of fields
        lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),1,j); hold on;
        h=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % Plot
            try
            shadedErrorBar(1:length(AllDays_PlotLearning.(datamatfield).(syl).meanFF), AllDays_PlotLearning.(datamatfield).(syl).meanFF, ...
                AllDays_PlotLearning.(datamatfield).(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1)
            
            h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
            catch err
            end
        end
        % annotate
        try
        legend(h,FieldsList)
        catch err
        end
        title(['Mean Pitch, ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (hz)','FontSize',12,'FontWeight','bold')
        xlabel('days','FontSize',12,'FontWeight','bold')
        
        % WN lines
        Fn_AnnotateWNLines(plotWNdays,ylim)
    end
    if plotLMANmusc==1;
        lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd}]);
    end
    
    
    
    % PLOT pitch deviation and zscore - each with own axis.
    lt_figure; hold on;
    hfig1=[];
    hfig2=[];
    for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        % PLOT MEAN PITCH SUBTRACTING BASELINE
        lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
        h=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % compile pitch means and SEM for each day, for this syl.
            try
                shadedErrorBar(1:length(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase), AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase, ...
                AllDays_PlotLearning.(datamatfield).(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1);
            
            h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
            catch err
            end
        end
        try
        legend(h,FieldsList)
        catch err
        end
        title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
        xlabel('days','FontSize',12,'FontWeight','bold')
        
        lt_plot_zeroline
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
        
        % PLOT Z-score
        lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,j*2); hold on;
        h=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % Plot
            try               
            plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'-','Color',plot_colors{jj},'LineWidth',1.5)
            h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
            catch err
            end
        end
        
        % annotate
        try
        legend(h,FieldsList)
        catch err
        end
            title(['Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (Z-scored to baseline)','FontSize',12,'FontWeight','bold')
        xlabel('days','FontSize',12,'FontWeight','bold')
        
        lt_plot_zeroline;
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
    end
    if plotLMANmusc==1;
        lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd}]);
    end
    
    
    
    % PLOT pitch deviation and zscore - making yaxis similar
    lt_figure; hold on;
    hfig1=[];
    hfig2=[];
    for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
        FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        % PLOT MEAN PITCH SUBTRACTING BASELINE
        hfig1(j)=lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,-1+j*2); hold on;
        h=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % compile pitch means and SEM for each day, for this syl.
            try
            shadedErrorBar(1:length(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase), AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase, ...
                AllDays_PlotLearning.(datamatfield).(syl).semFF,{'Color',plot_colors{jj},'LineWidth',1.5},1);
            
            h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_DevFromBase,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
            catch err
            end
            
        end
        try
        legend(h,FieldsList)
        catch err
        end
        title(['Day mean Pitch (minus baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (hz) (rel to baseline)','FontSize',12,'FontWeight','bold')
        xlabel('days','FontSize',12,'FontWeight','bold')
        
        lt_plot_zeroline
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
        
        % PLOT Z-score
        hfig2(j)=lt_subplot(length(Params.SeqFilter.SylLists.FieldsToPlot),2,j*2); hold on;
        h=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            % Plot
            try
            plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'-','Color',plot_colors{jj},'LineWidth',1.5);
            h(jj)=plot(AllDays_PlotLearning.(datamatfield).(syl).meanFF_zscore,'o','Color',plot_colors{jj},'MarkerFaceColor',plot_colors{jj},'MarkerSize',6);
            catch err
            end
            end
        
        % annotate
        try
        legend(h,FieldsList)
        catch err
        end
        title(['Z-scored mean pitch (relative to baseline), ' num2str(FirstDay) ' to ' num2str(LastDay)]);
        ylabel('FF (Z-scored to baseline)','FontSize',12,'FontWeight','bold')
        xlabel('days','FontSize',12,'FontWeight','bold')
        
        lt_plot_zeroline;
        Fn_AnnotateWNLines(plotWNdays,ylim)
        
        
        
    end
    if plotLMANmusc==1;
        lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd}]);
    end
    
    
    linkaxes(hfig1,'y');
    linkaxes(hfig2,'y');
    
end

%% PLOT SYLLABLES IN SEQUENCE ORDER - does both PBS and MUSC


if plotLMANmusc==1; % IF LMAN DATA EXISTS, THEN PLOT SEPARATE FOR MUSC AND PBS
    epochdata_fields={'EpochData', 'EpochData_MUSC'};
    datamatrix_fields={'DataMatrix', 'DataMatrix_MUSC'};
    
else
    epochdata_fields={'EpochData'};
    datamatrix_fields={'DataMatrix'};
end

for dd=1:length(epochdata_fields);
    epochdatafield=epochdata_fields{dd};
    datamatrixfield=datamatrix_fields{dd};
    
    % 1st and last WN bins
    try
    FirstWNBin=AllDays_PlotLearning.(epochdatafield).AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).FirstWNInd;
    LastWNBin=AllDays_PlotLearning.(epochdatafield).AllDaysSliding.(Params.PlotLearning.DayBinsFieldname).(syl).LastWNInd;
    catch err
        continue;
    end
    
    for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        lt_figure; hold on;
        
        X=length(Params.SeqFilter.SylLists.FieldsInOrder{i}); % how many syls?
        Y=AllDays_PlotLearning.(epochdatafield).MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline;
        
        Ymin=min(min(Y));
        Ymax=max(max(Y));
        Ylimits{i}=[Ymin-10 Ymax+10];
        Ylimits_zoom{i}=Ylimits{i}./5;
        
        for ii=1:X;
            lt_subplot(X,2,(-1+(X-ii+1)*2)); hold on;
            bar(Y(:,ii));
            try 
                ylim(Ylimits{i})
            catch err
            end
            ylabel(AllDays_PlotLearning.(epochdatafield).MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder{ii});
            
            % annotate WN bins
            line([FirstWNBin-0.5 FirstWNBin-0.5], ylim);
            line([LastWNBin-0.5 LastWNBin-0.5], ylim);
            
            
            if ii==X;
                xlabel('Day bin #');
            end
        end
        
        
        for ii=1:X;
            lt_subplot(X,2,((X-ii+1)*2)); hold on;
            bar(Y(:,ii));
            try 
                ylim(Ylimits_zoom{i})
            catch err
            end
            
            ylabel(AllDays_PlotLearning.(epochdatafield).MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder{ii});
            if ii==X;
                xlabel('Day bin #');
            end
            
            % annotate WN bins
            line([FirstWNBin-0.5 FirstWNBin-0.5], ylim);
            line([LastWNBin-0.5 LastWNBin-0.5], ylim);
            
        end
        
        if plotLMANmusc==1;
            lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd} ': Mean in 3-day bin']);
        else
            lt_subtitle('Mean FF in 3-day bin');
        end
    end
    
    % PLOT EACH DAY INDIVIDUALLY
    
    for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{i};
        
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        lt_figure; hold on;
        for ii=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{ii}; % actual syl name (e.g. 'a')
            lt_subplot(length(FieldsList),2,-1+(length(FieldsList)-ii+1)*2); hold on;
            
            % Plot pitch means and SEM for each day, for this syl.
            bar(AllDays_PlotLearning.(datamatrixfield).(syl).meanFF_DevFromBase);
            
            try
            ylim(Ylimits{i})
            catch err
            end
            
            ylabel(syl);
            if ii==length(FieldsList);
                xlabel('3-day bin #');
            end
            
            % WN lines
            Fn_AnnotateWNLines(plotWNdays,ylim)
        end
        
        % ZOOM IN
        for ii=1:length(FieldsList);
            syl=FieldsList{ii}; % actual syl name (e.g. 'a')
            lt_subplot(length(FieldsList),2,(length(FieldsList)-ii+1)*2); hold on;
            bar(AllDays_PlotLearning.(datamatrixfield).(syl).meanFF_DevFromBase);
            
            try
            ylim(Ylimits_zoom{i})
            catch err
            end
            
            ylabel(syl);
            if ii==length(FieldsList);
                xlabel('3-day bin #');
            end
            % WN lines
            Fn_AnnotateWNLines(plotWNdays,ylim)
        end
        
        if plotLMANmusc==1;
            lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd} ': (Syls in order) Mean learning per day']);
        else
            lt_subtitle('Syls in order) Mean learning per day');
        end
    end
end



%% PLOT HEAT MAP - currently only does PBS

% LEARNING (Hz from baseline) - 3 day bins

try 
figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(Params.PlotLearning.DayBinSize)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    % WN Lines
    line(xlim,[FirstWNBin-0.5 FirstWNBin-0.5])
    line(xlim,[LastWNBin-0.5 LastWNBin-0.5])

    %     global title
    subtitle('Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end

% ZOOM IN
figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Binned days during WN, Binsize: ' num2str(Params.PlotLearning.DayBinSize)]);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    try
    caxis(Ylimits_zoom{i})
    catch err
    end
    
    
    line(xlim,[FirstWNBin-0.5 FirstWNBin-0.5])
    line(xlim,[LastWNBin-0.5 LastWNBin-0.5])

    %     global title
    subtitle('ZOOM: Learning (mean Hz minus baseline) for each syllable - 3 day bins.');
end





% PLOT SAME, but day by day
figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    % WN start and end line
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('Learning (mean Hz minus baseline) for each syllable.');
end

% ZOOM IN

figure; hold on;
for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
    subplot(1,length(Params.SeqFilter.SylLists.FieldsInOrder),i); hold on;
    
    % UNDIR
    imagesc(DataMatrix.MatrixOverDaysforSylLists.FieldsInOrder{i}.meanFF_minusBaseline);
    
    colormap('gray');
    hbar=colorbar;
    ylabel(['Days']);
    Xtick=1:length(EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder); % one tick for each syl. needed.
    set(gca,'XTick',Xtick);
    set(gca,'XTickLabel',EpochData.MatrixOverDaysforSylLists.FieldsInOrder{i}.FieldNamesInOrder)
    
    try
    caxis(Ylimits_zoom{i})
    catch err
    end
    
        
    % WN start and end line
    
    line(xlim,[WNTimeOnInd-0.5 WNTimeOnInd-0.5])
    line(xlim,[WNTimeOffInd-0.5 WNTimeOffInd-0.5])
    
    %     global title
    subtitle('ZOOM: Learning (mean Hz minus baseline) for each syllable.');
end


catch err
end


%% PLOT MEAN OF LAST FEW DAYS IN ORDER OF SYLLABLES (and snapshot days, if desired) - currently only for PBS

try
% PLOT LEARNING AS MINUS BASELINE (HZ);
binfield=Params.PlotLearning.DayBinsFieldname;
lastday=EpochData.AllDaysSliding.(binfield).(syl).LastWNInd;

% only perform if there is data for last day
if length(EpochData.AllDaysSliding.(binfield).(syl).meanFF_minusBaseline)>=lastday & ...
        lastday>0;
    
    Y={};
    Ysem={};
    for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        % Get all syls
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
        plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        
        Y{ll}=[];
        Ysem{ll}=[];
        for jj=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            
            Y{ll}=[Y{ll} EpochData.AllDaysSliding.(binfield).(syl).meanFF_minusBaseline(lastday)];
            Ysem{ll}=[Ysem{ll} EpochData.AllDaysSliding.(binfield).(syl).semFF(lastday)];
        end
        
        Ymax(ll)=max(Y{ll});
        Ymin(ll)=min(Y{ll});
        
    end
    
    
    % Plot
    lt_figure; hold on;
    
    for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
        lt_subplot(length(Params.SeqFilter.SylLists.FieldsInOrder),1,ll); hold on;
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
        
        bar(Y{ll})
        
        try
        ylim([min(Ymin)-20 max(Ymax)+20]);
        catch err
        end
        
        % annotate
        Xtick=1:length(FieldsList); % one tick for each syl. needed.
        set(gca,'XTick',Xtick);
        set(gca,'XTickLabel',FieldsList,'FontSize',12,'FontWeight','bold')
        ylabel('FF (hz), minus baseline','FontSize',12,'FontWeight','bold')
        
    end
    lt_subtitle(['Learning (Hz minus baseline) last 3 WN days']);
else
    disp('Note, not plotting summary of last WN day because lacks data!!');
end
catch err
end
    
%% TO DO: CONVERT TO FUNCTION TO PLOT ANY DESIRED DAYS
% IF WANT TO PLOT A SNAPSHOT - need to designate in params
% Params.SeqFilter.DaysForSnapshot{1}={'09Dec2014','11Dec2014'};
% Params.SeqFilter.DaysToMark= {'11Dec2014-2400'}; % will mark all plots with lines here;

% only suppports one entry for days ATM.

if isfield(Params.SeqFilter,'DaysForSnapshot');
    Params.PlotLearning.fnamedays=(['from' Params.SeqFilter.DaysForSnapshot{1}{1} 'to' Params.SeqFilter.DaysForSnapshot{1}{2}]); % field name for structure below
    
    [A]=lt_convert_EventTimes_to_RelTimes(FirstDay,Params.SeqFilter.DaysForSnapshot{1});
    
    fday=A.JustDays_rel(1);
    lday=A.JustDays_rel(2);
    
    for i=1:length(SylFields_Unique);
        syl=SylFields_Unique{i};
        
        X=[];
        for ii=fday:lday; % over the days
            if ~isempty(AllDays_RawDatStruct{ii});
                if ~isfield(AllDays_RawDatStruct{ii}.data, syl);
                   % somtimes no data for this syl for this day...
                   continue
                end
            X=[X; cell2mat(AllDays_RawDatStruct{ii}.data.(syl)(:,1))];
            else % if empty, then tell user
                disp(['day ' num2str(ii) ' lacks data, but is trying to be used in Snapshot']);
            end
        end
        
        Snapshot.(Params.PlotLearning.fnamedays).(syl).FFvals=X;
        Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF=mean(X);
        Snapshot.(Params.PlotLearning.fnamedays).(syl).stdFF=std(X);
        Snapshot.(Params.PlotLearning.fnamedays).(syl).n=length(X);
        Snapshot.(Params.PlotLearning.fnamedays).(syl).sem=std(X)/sqrt(length(X)-1);
        
        % DEVIATION FROM BASELINE
        Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline=...
            Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF-EpochData.Baseline.(syl).meanFF;
        
        %z-score
        Snapshot.(Params.PlotLearning.fnamedays).(syl).Zscore=...
            Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline/EpochData.Baseline.(syl).stdFF;
        
        %             % GENERALIZATION STATS
        %         for i=1:length(SylFieldsAll);
        %             syl=SylFieldsAll{i};
        %
        %             for k=1:length(SylLists.TargetSyls);
        %                 targsyl=SylLists.TargetSyls{k}; % target;
        %
        %                 EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingHz=...
        %                     EpochData.Snapshot.(fnamedays).UNDIR.(syl).meanFF_minusBaseline/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).meanFF_minusBaseline;
        %
        %                 EpochData.Snapshot.(fnamedays).UNDIR.(syl).GeneralizationFrom.(targsyl).UsingZ=...
        %                     EpochData.Snapshot.(fnamedays).UNDIR.(syl).Zscore/EpochData.Snapshot.(fnamedays).UNDIR.(targsyl).Zscore;
        %             end
        %         end
        
    end

        % NOW PLOT
        % First, gather data in plotting format for BAR
        Y={};
        Ysem={};
        for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
            % Get all syls
            FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
            plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
            
            Y{ll}=[];
            Ysem{ll}=[];
            for jj=1:length(FieldsList); % how many fields within this set?
                syl=FieldsList{jj}; % actual syl name (e.g. 'a')
                
                Y{ll}=[Y{ll} Snapshot.(Params.PlotLearning.fnamedays).(syl).meanFF_minusBaseline];
                Ysem{ll}=[Ysem{ll} Snapshot.(Params.PlotLearning.fnamedays).(syl).sem];
            end
            
            Ymax(ll)=max(Y{ll});
            Ymin(ll)=min(Y{ll});
            
        end
        
        
        % Second, Plot
        figure; hold on;
        
        for ll=1:length(Params.SeqFilter.SylLists.FieldsInOrder);
            subplot(length(Params.SeqFilter.SylLists.FieldsInOrder),1,ll); hold on;
            FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{ll};
            
            bar(Y{ll})
            
            ylim([min(Ymin)-20 max(Ymax)+20]);
            
            % annotate
            Xtick=1:length(FieldsList); % one tick for each syl. needed.
            set(gca,'XTick',Xtick);
            set(gca,'XTickLabel',FieldsList,'FontSize',12,'FontWeight','bold')
            ylabel('FF (hz), minus baseline','FontSize',12,'FontWeight','bold')
            
        end
        [~,h]=subtitle(['Learning (Hz minus baseline) from ' Params.SeqFilter.DaysForSnapshot{1}{1} ' to ' Params.SeqFilter.DaysForSnapshot{1}{2}]);
        set(h,'FontSize',12,'FontWeight','bold')


AllDays_PlotLearning.EpochData.Snapshot=Snapshot;
end



%% PLOT PITCH CONTOUR OVER LEARNING

if plotLMANmusc==1; % IF LMAN DATA EXISTS, THEN PLOT SEPARATE FOR MUSC AND PBS
    epochdata_fields={'EpochData', 'EpochData_MUSC'};
    rawdata_fields={'data', 'data_MUSC'};
    
else
    epochdata_fields={'EpochData'};
    rawdata_fields={'data'};
end

for dd=1:length(epochdata_fields);
    epochdatafield=epochdata_fields{dd};
    rawdatafield=rawdata_fields{dd};
    
    DaysToPlot=[];
    % ---- Get subset of days to plot pitch contour of
    % last baseline day
    DaysToPlot(1)=Params.SeqFilter.BaselineDays(end);
    
    % 4 days into learning
    if Params.SeqFilter.BaselineDays(end)+4<NumDays;
    DaysToPlot(2)=Params.SeqFilter.BaselineDays(end)+4;
    else
    DaysToPlot(2)=NumDays;
    end
        
    % last day of learning
    DaysToPlot(3)=Params.PlotLearning.WNTimeOffInd;
    
    % any specified day to mark
    if isfield(Params.PlotLearning,'DaysToMarkInds');
        if ~isempty(Params.PlotLearning.DaysToMarkInds);
            DaysToPlot=[DaysToPlot cell2mat(Params.PlotLearning.DaysToMarkInds)];
        end
    end
    
    % reorder
    DaysToPlot=sort(DaysToPlot);
    % ---
    
    plotcolors=lt_make_plot_colors(length(DaysToPlot),1,[1 0.2 0.2]);
    
    % PLOT
    for j=1:length(Params.SeqFilter.SylLists.FieldsInOrder); % how many sets of fields (i.e. syls)?
        FieldsList=Params.SeqFilter.SylLists.FieldsInOrder{j};
        
        %     plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
        lt_figure; hold on;
        
        for jj=1:length(FieldsList); % how many fields within this set?
            
            syl=FieldsList{jj}; % actual syl name (e.g. 'a')
            subplot(1,length(FieldsList),jj); hold on;
            title(['Syl: ' syl]);
            
            % extract baseline PC mean and std
            baselinePC_mean=AllDays_PlotLearning.(epochdatafield).Baseline.(syl).pitchcontour_mean;
            baselinePC_std=AllDays_PlotLearning.(epochdatafield).Baseline.(syl).pitchcontour_STD;
            
            
            for jjj=1:length(DaysToPlot);
                dayInd=DaysToPlot(jjj);
                
                if length(AllDays_RawDatStruct)<dayInd;
                    continue; 
                end
                    
                if ~isempty(AllDays_RawDatStruct{dayInd})
                    if ~isfield(AllDays_RawDatStruct{dayInd}.(rawdatafield), syl);
                        continue;
                    end
                    if ~isempty(AllDays_RawDatStruct{dayInd}.(rawdatafield).(syl));
                        
                        % === calcualte deviation from baseline PC
                        % -- mean PC
                        PCmat=cell2mat(AllDays_RawDatStruct{dayInd}.(rawdatafield).(syl)(:,2));
                        PC_mean=mean(PCmat,1);
                        
                        % make sure the PCs are same length. if not, assume
                        % the longer one should be clipped at the end (i.e.
                        % the beginnings are the same)
                        if length(PC_mean)>length(baselinePC_mean);
                            PC_mean=PC_mean(1:length(baselinePC_mean));
                        end
                        
                        PC_minusBase=PC_mean-baselinePC_mean;
                        PC_zscore=PC_minusBase./baselinePC_std;
                        
                        % == get time windows
                        % what was time window used to get FF?
                        % -- use Recalculated FF if they exist:
                        if isfield(Params, 'RecalculateFF');
                            time_window=Params.RecalculateFF.pc_time_window_list(:,strcmp(Params.PlotLearning.SylFieldsAll, syl));
                        else
                            % use the original time window
                            time_window=Params.SeqFilter.pc_time_window_list{dayInd}(:,strcmp(fieldnames(AllDays_RawDatStruct{dayInd}.data), syl));
                        end
                        
                        
                        % get +/- 20ms
                        tbin_size=Params.DayRawDat.pc_T{1}(3)-Params.DayRawDat.pc_T{1}(2); % s in one time bin
                        
                        time_window_mod(1)=time_window(1)-floor(0.02/(tbin_size));
                        if time_window_mod(1)<1;
                            time_window_mod(1)=1;
                        end
                        
                        time_window_mod(2)=time_window(2)+floor(0.02/(tbin_size));
                        
                        % -- Plot
                        
                        try
                            % mean PC
                            % NOTE: do not use summary stats - instead, recalcualte
                            % mean as in lines above.  reason: muscimol data might
                            % be mixed with PBS.
                            %                     plot(AllDays_RawDatStruct{dayInd}.summary_stats.(syl).pitchcountour_mean,'Color',plotcolors{jjj},'LineWidth',2);
                            %                     text(1,AllDays_RawDatStruct{dayInd}.summary_stats.(syl).pitchcountour_mean(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                            
                            % deviation of PC from baseline
                            %                     plot(PC_minusBase,'Color',plotcolors{jjj},'LineWidth',2);
                            %                     text(1,PC_minusBase(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                            
                            % zscore of PC
                            plot(Params.DayRawDat.pc_T{1}(time_window_mod(1):time_window_mod(2)), PC_zscore(time_window_mod(1):time_window_mod(2)),'Color',plotcolors{jjj},'LineWidth',2);
                            
                            text(1,PC_zscore(end),['Day: ' num2str(dayInd)],'Color',plotcolors{jjj});
                            
                            % line indicating where got FF
                            line([Params.DayRawDat.pc_T{1}(time_window(1)) Params.DayRawDat.pc_T{1}(time_window(1))], ylim);
                            line([Params.DayRawDat.pc_T{1}(time_window(2)) Params.DayRawDat.pc_T{1}(time_window(2))], ylim);
                            
                            
                        catch err
                            disp('problem with pitch contour, had to use try');
                        end
                    end
                    
                end
            end
            
        end
        if plotLMANmusc==1;
            lt_subtitle([Params.PlotLearning.ExptCondition_codes{dd} ': Day mean pitch contour']);
        else
            lt_subtitle(['Day mean pitch contour']);
        end
    end
end


%% PLOT HIT RATE OVER LEARNING

% collect hit statuses for all trials
lt_figure; hold on;
title('Hit Rate each day (likely for catch songs) for all syls');
xlabel('day');
ylabel('Fraction of renditions hit');

plotcols=lt_make_plot_colors(length(SylFields_Unique),0,0);

for j=1:length(SylFields_Unique);
    syl=SylFields_Unique{j};
    
    hitrate=nan(1,NumDays); % to collect values to plot
    for i=1:NumDays;
        
        if ~isempty(AllDays_RawDatStruct{i});
            
            if ~isfield(AllDays_RawDatStruct{i}.data, syl);
                continue;
            end
            % collect data
            AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{i}=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,9));
            
            hitrate(i)=sum(AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{i})/length(AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{i});
        end
        
        % plot
        tmp=length(syl);
        
        if tmp==1;
            hplot(j)=plot(hitrate,'--','Color',plotcols{j},'MarkerFaceColor',plotcols{j});
        elseif tmp>1
            hplot(j)=plot(hitrate,'-o','Color',plotcols{j},'MarkerFaceColor',plotcols{j},'LineWidth',2);
        end
        
    end
end
legend(hplot,SylFields_Unique);


Fn_AnnotateWNLines(plotWNdays,ylim)


%% IS THERE CIRCADIAN FLUCTUATION THAT SHOULD BE SUBTRACTED FROM DATA?

lt_figure; hold on;
hplot=[];
hsplot=[];
syllist=SylFields_Unique;

plotcols=lt_make_plot_colors(length(syllist),0,0);

for j=1:length(syllist);
    syl=syllist{j};
    
    hsplot(j)=lt_subplot(ceil((length(syllist)/3)),3,j); hold on;
    title(syl);
    
    % collect values    
    FFvals_tot=[];
    Tvals_tot=[];
    for i=Params.SeqFilter.BaselineDays;
        if ~isempty(AllDays_RawDatStruct{i});
            if ~isfield(AllDays_RawDatStruct{i}.data, syl);
                continue;
            end
            FFvals=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,1));
            tmp=cell2mat(AllDays_RawDatStruct{i}.data.(syl)(:,6));
            [~, tmp] =lt_convert_datenum_to_hour(tmp);
            Tvals=tmp.hours;
            
            % collect vals across days
            FFvals_tot=[FFvals_tot; FFvals];
            Tvals_tot=[Tvals_tot; Tvals];
        end
    end
    
    % plot
    if isempty(FFvals_tot);
        continue;
    end
      
    hplot(j)=plot(Tvals_tot,FFvals_tot,'.','Color',plotcols{j});
    
    % get running avg
    [Tvals_tot, ind]= sort(Tvals_tot);
    FFvals_tot=FFvals_tot(ind);
    
    Tvals_tot_sm=lt_running_stats(Tvals_tot,10);
    FFvals_tot_sm=lt_running_stats(FFvals_tot,10);
    
    plot(Tvals_tot_sm.Mean,FFvals_tot_sm.Mean,'-k','LineWidth',2);
    
    xlim([6 22]);
end

legend(hplot,syllist);
linkaxes(hsplot,'x');
lt_subtitle('Baseline circadian fluctuation of pitch');


%% PLOT TRIAL BY TRIAL
if (0) % SKIPPING BECAUSE TAKES  ALONG TIME.

% TO DO - concatenate multiple days
% - overlay targ syl

SmoothBin=10; % 10 renditions
targsyl=Params.SeqFilter.SylLists.TargetSyls{1};

% PLOT pitch deviation and zscore - making yaxis similar
for j=1:length(Params.SeqFilter.SylLists.FieldsToPlot); % how many sets of fields (i.e. syls)?
    FieldsList=Params.SeqFilter.SylLists.FieldsToPlot{j};
    
    plot_colors=lt_make_plot_colors(length(FieldsList),0); % initiate colors for plot
    
    hfig1=lt_figure; hold on % for raw values with times
    hfig2=lt_figure; hold on; % for values by rendition
    hfig3=lt_figure; hold on; % for binned (early day vs. late);
    
    hsplot1=[];
    hsplot2=[];
    hsplot3=[];
    
    for jj=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{jj}; % actual syl name (e.g. 'a')
        
        % ================= combine raw data across all days
        FFvals_alldays=[];
        Tvals_alldays=[];
        FFvals_HalfBinned_alldays_mean=nan(NumDays,2); % col1=half1, col2=half2, rows = days
        FFvals_HalfBinned_alldays_sem=nan(NumDays,2); % col1=half1, col2=half2, rows = days
        
        for jjj=1:NumDays;
            dayInd=jjj;
            
            % skip if this syl has no data for today
            
            if isempty(AllDays_RawDatStruct{dayInd});
                continue;
            end
            
            if ~isfield(AllDays_RawDatStruct{dayInd}.data, syl);
                continue;
            end
            
            if isempty(AllDays_RawDatStruct{dayInd}.data.(syl));
                continue;
            end
            
            
            % ===== COLLECT FFVALS AND TVALS FOR TODAY
            % ffvals
            FFvals=cell2mat(AllDays_RawDatStruct{dayInd}.data.(syl)(:,1));
            
            % get times
            tmp=cell2mat(AllDays_RawDatStruct{dayInd}.data.(syl)(:,6));
            tmp=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,tmp);
            Tvals=tmp.FinalValue; % get times in form of (1.5 = day 1, noon)
            
            % === PUT INTO ACROSS DAYS ARRAY
            FFvals_alldays=[FFvals_alldays; FFvals];
            Tvals_alldays=[Tvals_alldays; Tvals];
            
            % === PUT INTO DAY BINNED (HALF VS. HALF) ARRAY
            inds1=1:floor(length(FFvals)/2);
            inds2=floor(length(FFvals)/2)+1:length(FFvals);
            
            FFval_mean_early=mean(FFvals(inds1));
            FFval_sem_early=lt_sem(FFvals(inds1));
            
            FFval_mean_late=mean(FFvals(inds2));
            FFval_sem_late=lt_sem(FFvals(inds2));
            
            FFvals_HalfBinned_alldays_mean(jjj,:) = [FFval_mean_early FFval_mean_late];
            FFvals_HalfBinned_alldays_sem(jjj,:)= [FFval_sem_early FFval_sem_late];
        end
        
        
        % ==== subtract (baseline) from FF values
        if ~isfield(AllDays_PlotLearning.EpochData.Baseline, syl);
        continue; 
        end
       
        baselineFF=AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
        %                 tmp=round(length(FFvals_alldays)/6);
        %                 baselineFF=mean(FFvals_alldays(1:tmp));
        FFvals_alldays=FFvals_alldays-baselineFF;
        % smooth
        FFvals_sm=lt_running_stats(FFvals_alldays,SmoothBin);
        Tvals_sm=lt_running_stats(Tvals_alldays,SmoothBin);
        
        % ========================================== Plot
        % === PLOT 1 - raw values with real times
        figure(hfig1);
        hsplot1(jj)=lt_subplot(ceil(length(FieldsList)/2),2,jj); hold on;
        title(syl);
        
        % RAW VALUES
        lt_plot(Tvals_alldays, FFvals_alldays, {'Color', plot_colors{jj}});
        % SMOOTHED
        if length(Tvals_sm.Mean)>1;
        shadedErrorBar(Tvals_sm.Mean,FFvals_sm.Mean,FFvals_sm.SEM,{'-','Color','k','LineWidth',6},1);
        %                 plot(Tvals_sm.Mean,FFvals_sm.Mean,'.','Color',plot_colors{jj})
        end
        
        % plot line for day mean
        plot(xlim,[0 0],'--');
        Xlim=xlim;
        text(Xlim(1),0,num2str(baselineFF));
        
        % Plot vertical lines
        line([Params.SeqFilter.BaselineDays(end)+1 Params.SeqFilter.BaselineDays(end)+1], ylim);
        line([Params.PlotLearning.WNTimeOffInd+1 Params.PlotLearning.WNTimeOffInd+1], ylim);
        
        
        % === PLOT 2 - by renditions
        figure(hfig2);
        hsplot2(jj)=lt_subplot(ceil(length(FieldsList)/2),2,jj); hold on;
        title(syl);
        
        X=1:length(FFvals_alldays);
        X_sm=lt_running_stats(X,SmoothBin);
        
        % RAW VALUES
        lt_plot(X, FFvals_alldays, {'Color', plot_colors{jj}});
        % SMOOTHED
                if length(X_sm.Mean)>1;
        shadedErrorBar(X_sm.Mean,FFvals_sm.Mean,FFvals_sm.SEM,{'-','Color','k','LineWidth',6},1);
                end
                
        % plot line for day mean
        plot(xlim,[0 0],'--');
        Xlim=xlim;
        text(Xlim(1),0,num2str(baselineFF));
        
        % --- Plot vertical lines
        % a line for each day
        for i=2:NumDays;
            % find the datapoint where the day transition is at
            % (i.e. the start of day i)
            tmp=Tvals_alldays-i;
            Ind=find(tmp>0, 1, 'first');
            if isempty(Ind)
                continue;
            end
            line([Ind-0.5 Ind-0.5],ylim);
            Ylim=ylim;
            text(Ind-10, 0, num2str(i-1), 'FontSize' ,13,'FontWeight','bold');
        end
        
        % ======== PLOT 3 - bin each day into first and second half
        figure(hfig3);
        hsplot3(jj)=lt_subplot(ceil(length(FieldsList)/2),2,jj); hold on;
        title(syl);
        
        % plot each day as 2 dots
        for n=1:NumDays;
            X=[n-0.2 n+0.2];
            Ymean=FFvals_HalfBinned_alldays_mean(n,:);
            Ysem=FFvals_HalfBinned_alldays_sem(n,:);
            
            errorbar(X, Ymean, Ysem, 'o-','Color', plot_colors{jj},'LineWidth',3);
        end
        
        % plot line for baseline mean
        line(xlim, [baselineFF baselineFF], 'Color','k');
        
        Fn_AnnotateWNLines(plotWNdays,ylim);
    end
    
    %          lt_subtitle('(subtract 1st 6th of day data');
    linkaxes(hsplot1,'xy');
    linkaxes(hsplot2,'xy');
     linkaxes(hsplot3,'x');
   
end

end    

%% TO DO:
% at: 
% Compile data into matrix form, for all syllable lists
% (apart from fields in order and fields of interest above

if (0)
UndirOrDirList={'UNDIR','DIR'};
lt_compile_seq_dep_pitch_data_PLOTDirUndir_ComplSylLists
end




% %% LABEL ALL RAW DATA AS EITHER MUSCIMOL OR NOT (if necessary)
% 
% if plotLMANmusc==1;
%     
%     [Params, AllDays_RawDatStruct, DataMatrix]= lt_seq_dep_pitch_PlotLearning_Musc(Params, AllDays_RawDatStruct, DataMatrix);
% 
% end
% 

%% OUTPUT
% 
% AllDays_PlotLearning.EpochData=EpochData;
% AllDays_PlotLearning.DataMatrix=DataMatrix;
% AllDays_PlotLearning.DataMatrix_Targ=DataMatrix_Targ;
%% SAVE
if saveON==1;
    
    cd(Params.PlotLearning.savedir);
    
    save('Params','Params');
    save('AllDays_PlotLearning','AllDays_PlotLearning');
    if save_AllDays_RawDat==1;
        save('AllDays_RawDatStruct','AllDays_RawDatStruct','-v7.3');
    end
    
    % write a text file that tells you when files were made
    fid1=fopen(['DONE_PlotLearning_' timestampSv '.txt'],'w');
    fclose(fid1);
    
    try
        cd FIGURES/PlotLearning
    catch err
        mkdir FIGURES/PlotLearning;
        cd FIGURES/PlotLearning;
    end
    
    lt_save_all_figs
    
    cd ../../
    
end


end

%% VARIOUS SUNFUNCTIONS

function Fn_AnnotateWNLines(plotWNdays,ylim)

global WNTimeOnInd
global WNTimeOffInd
global DaysToMarkInds

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd+0.5 WNTimeOffInd+0.5],ylim,'LineStyle','--','Color','r')
    
    for i=1:length(DaysToMarkInds);
        line([DaysToMarkInds{i} DaysToMarkInds{i}],ylim,'LineStyle','--','Color','k');
    end
end
end
