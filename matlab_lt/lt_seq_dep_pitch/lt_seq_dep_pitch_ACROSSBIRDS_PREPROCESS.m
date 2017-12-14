function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PREPROCESS(ExperimentList, PARAMS)

%% LT 8/13/15 - Preprocessing for across birds analysis

%% EXTRACT EXPERIMENT DATA STRUCTURES
NumExpt=length(ExperimentList);
SeqDepPitch_AcrossBirds=struct('birds',[]); % holding structure

for i=1:length(ExperimentList);
    
    % skip if is empty
    if isempty(ExperimentList{i});
        continue;
    end
       
    
    birdname=ExperimentList{i}{1};
    ExptID=ExperimentList{i}{2};
    NumTargs=ExperimentList{i}{3};
    PlotLearningDir=ExperimentList{i}{4};
    ConsolStartDate=ExperimentList{i}{5};
    ConsolEndDate=ExperimentList{i}{6};
    RegularExpressionsList=ExperimentList{i}{7};
    LMANinactivated=ExperimentList{i}{9};
    
    % === extract num days thr increased (single context phase)
    NumDaysAdapticeThr=ExperimentList{i}{11};
    if isempty(NumDaysAdapticeThr)
    disp([birdname '-' ExptID '-  ' 'EMPTY']);
    else
    disp([birdname '-' ExptID '-  ' num2str(NumDaysAdapticeThr)]);
    end
    
%     if NumTargs>2;
%         sdafasdf
%     end
    
% == old method. now using numtargs annotation
%     if ~isempty(ExperimentList{i}{8});
%         MultiDir_OneDayBeforeStart=ExperimentList{i}{8}{1};
%         MultiDir_LastDay=ExperimentList{i}{8}{2};
%         
%         num_multidir_syls=length(ExperimentList{i}{8})-2;
%         MultiDirSyls={};
%         for j=1:num_multidir_syls;
%             MultiDirSyls{j}=ExperimentList{i}{8}{j+2};
%         end
%         
%         multidir_exists=1;
%     else
%         multidir_exists=0;
%     end
    

% +++++ EXTRACT DIFFERENT EPOCHS - don't care which order the epochs came
% in. just to compare across experiments
    
    % === Extract multidir,
        multidir_exists=0;
    if NumTargs==2 | NumTargs==0 | NumTargs==4
        MultiDir_OneDayBeforeStart=ExperimentList{i}{8}{1};
        MultiDir_LastDay=ExperimentList{i}{8}{2};
        
        num_multidir_syls=length(ExperimentList{i}{8})-2;
        MultiDirSyls={};
        for j=1:num_multidir_syls;
            MultiDirSyls{j}=ExperimentList{i}{8}{j+2};
        end
        
        multidir_exists=1;
    end
    
    if NumTargs==5
        MultiDir_OneDayBeforeStart=ExperimentList{i}{10}{1};
        MultiDir_LastDay=ExperimentList{i}{10}{2};
        
        num_multidir_syls=length(ExperimentList{i}{10}{3});
        MultiDirSyls={};
        for j=1:num_multidir_syls;
            MultiDirSyls{j}=ExperimentList{i}{10}{3}{j};
        end
        multidir_exists=1;
    end
    
    
    % ======= Extract SAME DIR    
    samedir_exists=0;
    if NumTargs==3 | NumTargs==5
        SameDir_OneDayBeforeStart=ExperimentList{i}{8}{1};
        SameDir_LastDay=ExperimentList{i}{8}{2};
        
        num_samedir_syls=length(ExperimentList{i}{8})-2;
        SameDirSyls={};
        for j=1:num_samedir_syls;
            SameDirSyls{j}=ExperimentList{i}{8}{j+2};
        end
        
        samedir_exists=1;
    end
    
    if NumTargs==4
        SameDir_OneDayBeforeStart=ExperimentList{i}{10}{1};
        SameDir_LastDay=ExperimentList{i}{10}{2};
        num_samedir_syls=length(ExperimentList{i}{10}{3});
        SameDirSyls={};
        for j=1:num_samedir_syls;
            SameDirSyls{j}=ExperimentList{i}{10}{3}{j};
        end
        samedir_exists=1;
    end
        
    
    
        
    
    
    % Slot into structure
    if isempty(SeqDepPitch_AcrossBirds.birds);
        % then no birds slotted yet
        % start bird 1
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.PlotLearningDir=PlotLearningDir;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.NumTargs=NumTargs;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.ExptID=ExptID;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.RegularExpressionsList=RegularExpressionsList;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.LMANinactivated=LMANinactivated;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.NumDaysAdapticeThr_col1OK_col2complete=NumDaysAdapticeThr;
        SeqDepPitch_AcrossBirds.birds{1}.birdname=birdname;
        
        
        if multidir_exists==1;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
        end
        
        if samedir_exists==1;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.SameDir_OneDayBeforeStart=SameDir_OneDayBeforeStart;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.SameDir_LastDay=SameDir_LastDay;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.SameDirSyls=SameDirSyls;
        end
        
        
        
    else
        % find the right bird to slot data into
        c=0;
        for j=1:length(SeqDepPitch_AcrossBirds.birds);
            if strcmp(SeqDepPitch_AcrossBirds.birds{j}.birdname,birdname);
                % then this is right bird
                % enter new data, new expt
                ind=length(SeqDepPitch_AcrossBirds.birds{j}.experiment)+1;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.PlotLearningDir=PlotLearningDir;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.NumTargs=NumTargs;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.ExptID=ExptID;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolStartDate=ConsolStartDate;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolEndDate=ConsolEndDate;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.RegularExpressionsList=RegularExpressionsList;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.LMANinactivated=LMANinactivated;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.NumDaysAdapticeThr_col1OK_col2complete=NumDaysAdapticeThr;
                SeqDepPitch_AcrossBirds.birds{j}.birdname=birdname;
                
                if multidir_exists==1;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_LastDay=MultiDir_LastDay;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.num_multidir_syls=num_multidir_syls;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.MultiDirSyls=MultiDirSyls;
                end
                
                if samedir_exists==1;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.SameDir_OneDayBeforeStart=SameDir_OneDayBeforeStart;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.SameDir_LastDay=SameDir_LastDay;
                    SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.SameDirSyls=SameDirSyls;
                end
                
                
                
                % change c to show that bird found.
                c=1;
                
                continue % quit loop
            end
        end
        if c==0
            % then no bird was found, start new bird
            ind=length(SeqDepPitch_AcrossBirds.birds)+1;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.PlotLearningDir=PlotLearningDir;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.NumTargs=NumTargs;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.ExptID=ExptID;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.RegularExpressionsList=RegularExpressionsList;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.LMANinactivated=LMANinactivated;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.NumDaysAdapticeThr_col1OK_col2complete=NumDaysAdapticeThr;
            SeqDepPitch_AcrossBirds.birds{ind}.birdname=birdname;
            
            if multidir_exists==1;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
            end
            
            if samedir_exists==1;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.SameDir_OneDayBeforeStart=SameDir_OneDayBeforeStart;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.SameDir_LastDay=SameDir_LastDay;
                SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.SameDirSyls=SameDirSyls;
            end
            
            
        end
    end
end



%-- LOAD AND COMPILE STRUCTURES
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        % load
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/Params']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params=tmp.Params;
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/AllDays_PlotLearning']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning=tmp.AllDays_PlotLearning;
        
    end
end


% -- HOW many birds? 
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


% -- MAKE RANDOM EDITS (ad hoc)
 
% pu37wh20, % remove bccB as this was actually also a target (along with dccB);
% Don't need to, not in data nymore...
% for i=1:NumBirds;
%     if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu37wh20');
%     
%         for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%             if strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
%                 
% %                 % remove stuff
% %                 ii=find(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsSame,'bccB'))
%                 
%             end
%     
%     end
% end



%% == MAKE SURE Unique syls are defined for all expts (and if not, then do so)

SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_UniqueSyls(SeqDepPitch_AcrossBirds);

%% Remove any syls (from list of unique syls)
% ===== THROW OUT BAD SYLS (remove from "Syls_Unique"); % will remove from
% both single dir and bidir analysis simultaneously.
% -- if expt does not exist, or syl does not exist, is fine.
close all;
if PARAMS.global.remove_bad_syls==1;
[SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_RemoveSyls(SeqDepPitch_AcrossBirds, PARAMS);
end


% ===== PLOT SUMMARY OF METADATA OF EXPERIMENTS (E.G. SAMPLE SIZE)
% IN PROGRESS
% e.g. bar plot, each col is a bird, with height as num experiments of
% different ytpes)


%% Correlations and Acoustic distance analyses.
close all;
plotPCA_stuff=PARAMS.preprocess.plotPCA_stuff;
Zscore_use_rends=PARAMS.preprocess.Zscore_use_rends;
load_old_vector_space=PARAMS.preprocess.load_old_vector_space; % default, 0, recreates vector space on each run. if 1, then loads previous version.
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_CORR_and_PCA(SeqDepPitch_AcrossBirds, PARAMS, plotPCA_stuff, Zscore_use_rends,load_old_vector_space);




%% PLOT AND EXTRACT LEARNING PARAMS, LEARNING, SUMMARY PLOTS FOR EACH EXPERIMENT
close all;
DayBinSize=PARAMS.global.ConsolBinDaySize;
[SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_SummaryPlot(SeqDepPitch_AcrossBirds, PARAMS, DayBinSize);


%% MAKE SURE ALL LMAN EXPERIMENTS HAVE LMAN RAW DATA (e.g. regexpr) - IF NOT, RUN ANALYSES

currdir=pwd;
close all;

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % SKIP IF THIS IS NOT LMAN EXPERIMENT
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0;
            continue;
        end
        
        
        % ==== FIRST, make sure that regexpr has baselkine data for both
        % ---- musc (must already have pbs, because checked in previous code)
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr, 'baseline_MUSC');
            % Then need to go and do those analysis
            
            close all;
            disp(' ')
            disp('NEED TO DO ANALYSIS REGEXPR for MUSC');
            disp('Running...')
            
            % 1) go to dir
            cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
            
            % 2) run RegExpr and RegExprPlot analyses
            clear Params;
            clear AllDays_RegExpr;
            clear AllDays_RawDatStruct;
            clear AllDays_PlotLearning
            
            load('Params');
            load('AllDays_RawDatStruct');
            load('AllDays_PlotLearning');
            
            % ------ RUN REG EXPR
            Params.RegExpr.expressions=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
            saveON=1;
            DoLMAN=1;
            [Params, AllDays_RegExpr] = lt_seq_dep_pitch_RegExpr(Params, AllDays_RawDatStruct, saveON, DoLMAN, AllDays_PlotLearning);
            
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            LMANon=1;
            SuppressQueries=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon, SuppressQueries);
            
            % ---------- SAVE OUTPUTS
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
        end
        
        % ==== SECOND, make sure has baseline data for PBS
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr, 'baseline');
            % Then need to go and do those analysis (does not remove any
            % anlsyis of LMAN data done above).
            
            close all;
            disp(' ')
            disp('NEED TO DO ANALYSIS REGEXPR for PBS');
            disp('Running...')
            
            % 1) go to dir
            cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
            
            % 2) run RegExpr and RegExprPlot analyses
            clear Params;
            clear AllDays_RegExpr;
            
            load('Params');
            load('AllDays_RegExpr');
            
            % ------ RUN REG EXPR PLOR
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            LMANon=0;
            SuppressQueries=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon, SuppressQueries);
            
            % ---------- SAVE OUTPUTS
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
        end
    end
end
disp('DONE!...');
cd(currdir);


%% FILTER - To give each syl a sylID (its features)
% close all; 
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_FILTER(SeqDepPitch_AcrossBirds, PARAMS);


%% [LMAN] FILTER - To give each syl a sylID (its features) (LMAN EXPERIMENTS)
% Modified to allow to extract for LMAN data too.  Only does one thing -
% takes the output Syl_ID_Dimensions field and adds that as a new field
% Syl_ID_Dimensions_LMAN to the structure

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexperiments;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            % Then this has LMAN data
            
            
            % === FIRST, extract only that experiment (don't overwrite
            % struct) (replace bird 1 expt 1 with the current expt)
            SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
            
            SeqDepPitch_AcrossBirds.birds{1}.experiment{1}=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii};
            
            filter = 'LMAN';
            [SeqDepPitch_AcrossBirds_LMAN, ~]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
            
            
            % ==== SECOND, change the name of the field so that the next
            % program analyzes LMAN data, instead of PBS. (also remove the
            % SylID field that was previosly done).
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}=rmfield(SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}, ...
                'Syl_ID_Dimensions');
            
            disp(' ');
            disp('NOTE: LMAN Data filter will work for: ');
            % ---- To get correlations
            disp('Correlations...');
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline=...
                SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC;
            
            disp('Acoustic distance');
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_StructureStats.data=...
                SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_StructureStats.data_MUSC;
            
            
            disp('DOES NOT YET WORK FOR: ');
            disp('Learning_by_target (not important)')
            disp(' ');
            
            % ==== THIRD, RUN FILTER on LMAN DATA
            [SeqDepPitch_AcrossBirds_LMAN, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_FILTER(SeqDepPitch_AcrossBirds_LMAN, PARAMS);
            
            
            
            % ===== FOURTH, Take the SylID field from filtered LMAN data and put it
            % intot he global all experimtns structure
            SeqDepPitch_AcrossBirds_ORIGINAL.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN=SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Syl_ID_Dimensions;
            
            % ====== RESTORE THE ORIGINAL SeqDepPitch_AcrossBirds, with one
            % change (Syl_ID_Dimensions_LMAN is added for one expt/bird);
            SeqDepPitch_AcrossBirds=SeqDepPitch_AcrossBirds_ORIGINAL;
        end
    end
end

%% BASELINE DATA - MAKE SURE ALL TVALS IN BASELINE DATA (I DID NOT CODE THIS INTO THE ORIGINAL PLOT LEARNING FUNCTION)

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        baseline_days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % ==== GET VECTOR OF BASELINE FF AND TIMES
            ffvals_all=[];
            tvals_all=[];
            
            for k=baseline_days;
                ffvals_all=[ffvals_all cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k})];
                tvals_all=[tvals_all cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{k})];
            end
            
            % ---- make sure is same as what was done before
            n_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).n;
            meanFF_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
            
            n_new=length(ffvals_all);
            meanFF_new=mean(ffvals_all);
            
            if n_old~=n_new || meanFF_old~=meanFF_new;
                disp('PROBLEM - NEW BASELINE ANALYSIS DOES NOT MATCH OLD BASELINE');
                asdfasdfasfasd;
            end
            
            % ----- UPDATE DATA
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals=tvals_all; % Tvals
            %             disp([num2str(i) ' ' num2str(ii) ' ' syl]);
            
            % ============ IF IS LMAN, ALSO GET WITHIN TIME WINDOW
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                ffvals_all=[];
                tvals_all=[];
                
                for k=baseline_days;
                    ffvals_all=[ffvals_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{k}];
                    tvals_all=[tvals_all SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{k}];
                end
                
                % ---- make sure is same as what was done before
                n_old=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                meanFF_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                
                n_new=length(ffvals_all);
                meanFF_new=mean(ffvals_all);
                
                if n_old~=n_new || meanFF_old~=meanFF_new;
                    disp('PROBLEM - NEW BASELINE ANALYSIS DOES NOT MATCH OLD BASELINE');
                    asdfasdfasfasd;
                end
                
                % ----- UPDATE DATA
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals_WithinTimeWindow=tvals_all; % Tvals
                %             disp([num2str(i) ' ' num2str(ii) ' ' syl]);
                
            end
            
            
        end
    end
end

disp('--')
disp('BASELINE DAT FFVALS AND TVALS CHECKED')


%% WHAT is direction of target learning?

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
            % ==== targ learning dir
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targ_learndir=sign(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).ConsolEnd_meanFF);

            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir=targ_learndir;
    end
end

%% ======== for LMAN experiments, extract zscore within time window
% this is needed for Zscore data, if I want to use Zscore data. but I am
% overwriting with last 2 baseline days (in other code), so don't need to do this.

for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            for j=1:length(SylsUnique)
                syl=SylsUnique{j};
                
                % -- baseline stats
                base_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                base_std=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                
                % -- calc day to day stats
                meanFFdev_alldays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow;
                zscore_alldays=meanFFdev_alldays./base_std;
                
                % == OUTPUT
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore_WithinTimeWindow=zscore_alldays;
                
              
            end
        end
        
    end
end


%% ZSCORE DATA AND DEFINE LEARNING (I'M USING ZSCORE)
% testing things
PARAMS.zscore.zthresh_learning=0.8; % take 1st N days after target passes this value of shift (will take max of default day (see below) or this day)
PARAMS.zscore.default_min_day=3; 
PARAMS.zscore.max_day=3; % inclusive, limit days to within min and max days (wn start =1) then throw out this experiment
PARAMS.zscore.N_days_post_thresh=2; % how many days to use to quantify learning
PARAMS.zscore.account_for_NoData=1; % If early WN days don't have data, then the window will be moved forward
% the number of no-data days (assumes no singing early days)
[SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_ZSCORE(SeqDepPitch_AcrossBirds, PARAMS);


%% DEFINE LEARNING - all future analysis will be done with this metric
% use z-score or FF.  which epoch?
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            
            % ============ DEFINE THE LEARNING METRIC
            if PARAMS.global.learning_metric=='zscore';
                % then use z-score, with epoch defined in ZSCORE function above
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh,'DATA'); % some don't have z-score defined becasue not enough elarning.
                    learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore;
                    learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).sem_Zscore;
                    learning_DayInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
                else
                    learning=nan;
                    learning_sem=nan;
                    learning_DayInds=nan;
                end
                
                
            elseif PARAMS.global.learning_metric=='ff_consolstart';
                % then use mean FFval at time of consolidation start
                
                learning = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.consolid_start_meanFF_minbase;
                learning_sem=lt_sem(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FFvals_minbase);
                consol_start=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
                consol_start_numdays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Consol_DayBins;
                learning_DayInds=consol_start:consol_start+consol_start_numdays-1;

            else
                disp('PROBLEM - LEARNING METRIC NOT DEFINED PROPERLY, ANALYSES WILL NOT WORK');
                pause;
            end
            
            % ---- PUT INTO OUTPUT STRUCTURE
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean=learning;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem=learning_sem;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.dayIndsUsed=learning_DayInds;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed=learning_DayInds;
            
        end
        
        
        % ======== GET REALTUIVE TO TARGSYL
        % --- targsyl info
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learning_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning_rel_targ=learning/learning_targ;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ=learning_rel_targ;
            
        end
    end
end


%% ==== make sure extracted hit rate infor for all LMAN experiments [not necessarily within time window]
% NOTE: will have to restart entire analysis if runnign this is required.

lt_seq_dep_pitch_ACROSSBIRDS_ExtrHitsMUSC(SeqDepPitch_AcrossBirds);
            
%% ===== HIT RATES - for LMAN experiments extract hit rates within time window [for both PBS and MUSC data]

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
DatFields={'DataMatrix_MUSC', 'DataMatrix'};

for datfield=DatFields;
    datfield=datfield{1};
    for i=1:NumBirds;
        numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexperiments;
            
            %  -- check if this expt is lman
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
                continue
            end
            
            syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            
            for j=1:length(syls_unique);
                syl=syls_unique{j};
                
                NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).FFvals);
                
                % ---- pad number of days with empty days if necessary
                if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).HitRateStatus)<NumDays;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).HitRateStatus{NumDays}=[];
                end
                    
                for day=1:NumDays
                    
                    hitRate_allsongs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).HitRateStatus{day};
                    tvals_allsongs=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).Tvals{day});
                    tvals_timewindow=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).Tvals_WithinTimeWindow{day};
                    
                    ffvals_allsongs=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).FFvals{day});
                    ffvals_timewindow=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).FFvals_WithinTimeWindow{day};
                    
                    indsToKeep = find(ismember(tvals_allsongs, tvals_timewindow));
                    
                    % --- debug
                    if length(indsToKeep) ~= length(ffvals_timewindow)
                        if strcmp(datfield, 'DataMatrix_MUSC')
                            % if only one is empty, then there is problem. likley
                            % because musc data is taken from songs that are labeled
                            % PBS. solution, check those songs
                            hitRate_allsongs=[hitRate_allsongs; SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{day}];
                            tvals_allsongs=[tvals_allsongs cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day})];
                            ffvals_allsongs=[ffvals_allsongs cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day})];
                            
                            indsToKeep = find(ismember(tvals_allsongs, tvals_timewindow));
                        else
                            disp('PROBLEM!! - trying to extract PBS data, what I previously did is likely to take some MUSC data and use as PBS, but that code not yet written');
                            cafaedweacacea3ra;
                        end
                    end
                    
                    if length(indsToKeep) ~= length(ffvals_timewindow)
                        disp('PROBLEM - not sure why inds to keep extracted is not same as ffvals timewindow');
                        afsceaea3rawvarvrav;
                    end
                    
                    if ~isempty(indsToKeep) & ~isempty(ffvals_timewindow)
                        assert(all(ffvals_allsongs(indsToKeep)==ffvals_timewindow), 'problem - mathing hit rates to songs might have problem');
                    end
                    
                    % ===== COLLECT HIT RATES
                    hitRate_timewindow=hitRate_allsongs(indsToKeep);
                    
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.(datfield).(syl).HitRateStatus_WithinTimeWindow{day}=hitRate_timewindow;
                    
                    
                    % -- debug
                    assert(length(hitRate_timewindow)==length(ffvals_timewindow), 'WHAT!!!???');
                    
                    
                end
                
            end
        end
    end
end

%% ==== all LMAN experiments, make sure baseline has MUSC data


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts
        
        % === skip if is not LMAN
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
            continue
        end
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            % === collect all baseline rends from day data
            BaselineDays=1:SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd-1;
            
            FFvalsMUSC=[];
            TvalsMUSC=[];
            for day=BaselineDays
                % --- if this day is musc, then take
                if ~isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day})
                    
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                    tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                    
                    FFvalsMUSC=[FFvalsMUSC ffvals];
                    TvalsMUSC=[TvalsMUSC tvals];
                    
                end
            end
            
            % ==== compare, make sure is same as baseline data that's
            % already there.
            % --- current baseline stuff
            baseRawFF=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow;
            assert(length(baseRawFF)==length(FFvalsMUSC), 'asdffasfsd');
            assert(mean(baseRawFF)==mean(FFvalsMUSC), 'asdffasfsd');
            assert(length(TvalsMUSC) == length(FFvalsMUSC), 'dasfasdf');
            

            % --- add tvals
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).Tvals_WithinTimeWindow=TvalsMUSC;
        end
        
        
        
        
    end
end
%% ==== make sure there are no empty corrs (replace with nan)
% 
% NumBirds=length(SeqDepPitch_AcrossBirds.birds);
% 
% for i=1:NumBirds
%     numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     for ii=1:numexpts
%         
%         SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         
%         
%         for j=1:length(SylsUnique)
%             
%             syl=SylsUnique{j};
%             
%             
%             
%             
%            
%         end
%         
%         
%     end
%     
%     
%     
% end
% 
%                     if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS, 'motif_by_motif');
%                         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, syl2);
%                             
%                             corrmotifPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
%                             corrmotifMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
%                         end
%                     end
%                     
%                     if isempty(corrmotifPBS)
%                         corrmotifPBS=nan;
%                     end
% 
% 


%% ==== EXTRACT LAST SINGLE DIR DAY
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ')
disp(' ========= LAST SINGLE DIR DAY ==== ');

for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        daybeforemulti=1000; % 1000 will always be larger than window
        daybeforeSame=1000;
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            daybeforemulti=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        end
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart_Ind')
            daybeforeSame=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind;
        end
        
        LastSingleDirDay=min([daybeforemulti daybeforeSame]);
        
        if LastSingleDirDay==1000;
            % then entire expt is single dir.
            % -- last day with data
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            lastDayWithDat=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals);
            
            % -- last WN day
            LastSingleDirDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            
            % -- take minimum
            LastSingleDirDay=min([LastSingleDirDay lastDayWithDat]);
        end
        
        
        
        % ===== output
        disp([birdname '-' exptname ': ' num2str(LastSingleDirDay)]);
        
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SingleDirLastDay=LastSingleDirDay;
        
    end
end


