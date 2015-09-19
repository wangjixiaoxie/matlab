function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PREPROCESS(ExperimentList, PARAMS);

%% LT 8/13/15 - Preprocessing for across birds analysis

%% EXTRACT EXPERIMENT DATA STRUCTURES
NumExpt=length(ExperimentList);
SeqDepPitch_AcrossBirds=struct('birds',[]); % holding structure

for i=1:length(ExperimentList);
    birdname=ExperimentList{i}{1};
    ExptID=ExperimentList{i}{2};
    NumTargs=ExperimentList{i}{3};
    PlotLearningDir=ExperimentList{i}{4};
    ConsolStartDate=ExperimentList{i}{5};
    ConsolEndDate=ExperimentList{i}{6};
    RegularExpressionsList=ExperimentList{i}{7};
    LMANinactivated=ExperimentList{i}{9};
    
    % Multidir stuff
    if ~isempty(ExperimentList{i}{8});
        MultiDir_OneDayBeforeStart=ExperimentList{i}{8}{1};
        MultiDir_LastDay=ExperimentList{i}{8}{2};
        
        num_multidir_syls=length(ExperimentList{i}{8})-2;
        MultiDirSyls={};
        for j=1:num_multidir_syls;
            MultiDirSyls{j}=ExperimentList{i}{8}{j+2};
        end
        
        multidir_exists=1;
    else
        multidir_exists=0;
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
        SeqDepPitch_AcrossBirds.birds{1}.birdname=birdname;
        
        
        if multidir_exists==1;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
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
                SeqDepPitch_AcrossBirds.birds{j}.birdname=birdname;
                
                        if multidir_exists==1;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.MultiDir_LastDay=MultiDir_LastDay;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.num_multidir_syls=num_multidir_syls;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.INFORMATION.MultiDirSyls=MultiDirSyls;
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
            SeqDepPitch_AcrossBirds.birds{ind}.birdname=birdname;
            
                    if multidir_exists==1;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_OneDayBeforeStart=MultiDir_OneDayBeforeStart;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.MultiDir_LastDay=MultiDir_LastDay;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.num_multidir_syls=num_multidir_syls;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.INFORMATION.MultiDirSyls=MultiDirSyls;
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


% CHECK DIRECTION OF PITCH SHIFT


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
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_CORR_and_PCA(SeqDepPitch_AcrossBirds, PARAMS);




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
        % musc (must already have pbs, because checked in previous code)
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
close all; 
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
            SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline=...
                SeqDepPitch_AcrossBirds_LMAN.birds{1}.experiment{1}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC;
            
            disp('Correlations...');
            
            disp('DOES NOT YET WORK FOR: ');
            disp('Acoustic distance');
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

%% BASELINE DATA - MAKE SURE ALL HAVE FFVALS AND TVALS IN BASELINE DATA (I DID NOT CODE THIS INTO THE ORIGINAL PLOT LEARNING FUNCTION)

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
            end
            
            % ----- UPDATE DATA 
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals=tvals_all; % Tvals
            disp([num2str(i) ' ' num2str(ii) ' ' syl]);
        end
    end
end
