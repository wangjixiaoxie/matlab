function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, ...
    PARAMS, RepeatsOnly, OnlyUseSylsInSylsUnique, DaysToPlot, ExcludeSeqLearning, ExcludeNotFullyLabeled,...
    ExcludeIfHasSameDirBeforeDiffDir, ExcludeIfFirstTargDriveMore, RawShowOnlySameType, UseOldVersion, ...
    ThrowOutIfNotEnoughData)
%% LT 12/27/15 - plots day to day shift at two targets. Are they constrained to shift together?
% ExcludeSeqLearning=1; % then skips pu53, since sequence learning
% ExcludeNotFullyLabeled = 1; ad hoc, remove experiemtns I notice don't
% have enough days. see code below
% ExcludeIfHasSameDirBeforeDiffDir = 1; 
if ~exist('ExcludeNotFullyLabeled', 'var');
    ExcludeNotFullyLabeled=0;
end

if ~exist('ExcludeIfHasSameDirBeforeDiffDir', 'var');
    ExcludeIfHasSameDirBeforeDiffDir=0;
end

if isempty('ThrowOutIfNotEnoughData')
    ThrowOutIfNotEnoughData = 0;
end

adHocPlotSameType = 1; % if 1, then plots same type, nontarg (i.e. context 3) each individual syl
% as opposed to the current code, which takes mean within each expt. If set
% to 0, then plots diff type.

%% EXTRACT ONLY EXPERIMENTS WITH MULTIDIR
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

expts_to_remove=[];
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    expts_to_remove{i}=[];
    
    for ii=1:NumExperiments
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
            expts_to_remove{i}=[expts_to_remove{i} ii];
        end
    end
end

% actually remove those experiments from structure
for i=1:length(expts_to_remove);
    if ~isempty(expts_to_remove{i});
        SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
        disp(['TO GET ONLY MULTIDIR EXPERIMENTS, REMOVED: ' num2str(i) '; expt: ' num2str(expts_to_remove{i})]);
    end
end


PARAMS.global.MULTIDIR.expts_to_remove=expts_to_remove;
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

% total num experiments
TotalNumExperiments=0;
for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexperiments;
        TotalNumExperiments=TotalNumExperiments+1;
    end
end



%% ==== KEEP ONLY EXPTS WHERE ALL MULTIDIR TARGS ARE PART OF UNIQUE SYLS.
if OnlyUseSylsInSylsUnique==1;
    NumBirds=length(SeqDepPitch_AcrossBirds.birds);
    
    expts_to_remove=[];
    for i=1:NumBirds;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        expts_to_remove{i}=[];
        for ii=1:NumExperiments
            
            MultiDirSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
            for j=1:length(MultiDirSyls);
                syl=MultiDirSyls{j};
                
                if ~any(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique, syl))
                    expts_to_remove{i}=[expts_to_remove{i} ii];
                    break
                end
                
                
            end
        end
    end
    
    % actually remove them
    disp('Removed because syl not in sylsunique:');
    for i=1:length(expts_to_remove);
        if ~isempty(expts_to_remove{i});
            
            for k=1:length(expts_to_remove{i});
                exptInd=expts_to_remove{i}(k);
                
                birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
                exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{exptInd}.ExptID;
                disp([birdname ' - ' exptname]);
            end
            SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
        end
    end
    
    
    
    PARAMS.global.MULTIDIR.expts_to_remove=expts_to_remove;
    NumBirds=length(SeqDepPitch_AcrossBirds.birds);
    
    % total num experiments
    TotalNumExperiments=0;
    for i=1:NumBirds;
        numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        for ii=1:numexperiments;
            TotalNumExperiments=TotalNumExperiments+1;
        end
    end
end


%% === REMOVE EXPTS WITH MORE THAN 2 TARGS

disp(' ------ ');
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

expts_to_remove=[];
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    expts_to_remove{i}=[];
    for ii=1:NumExperiments
        
        MultiDirSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        if length(MultiDirSyls)>2
            expts_to_remove{i}=[expts_to_remove{i} ii];
        end
    end
end

% actually remove them
disp('Removed because >2 multidir syls:');
for i=1:length(expts_to_remove);
    if ~isempty(expts_to_remove{i});
        
        for k=1:length(expts_to_remove{i});
            exptInd=expts_to_remove{i}(k);
            
            birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{exptInd}.ExptID;
            disp([birdname ' - ' exptname]);
        end
        SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
    end
end


%% Extract only repeats experiments (if desired)

if RepeatsOnly==1;
    
    filter = 'MultiDirAreInRepeats';
    SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
    
end

NumBirds=length(SeqDepPitch_AcrossBirds.birds);




%% COLLECT ALL DAY TO DAY TRAJECTORIES - USING DAYS LOCKED TO BIDIR DAY START.
SameType_Others = []; % ad hoc, for revisions
ExptcountAll = [];

DATSTRUCT=struct;
Exptcount=1;
for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        if ExcludeSeqLearning==1;
           if strcmp(birdname, 'pu53wh88') && strcmp(exptname, 'SeqDepPitchShift3');
               disp(['SKIPPING ' birdname '- ' exptname ': sequence learning']);
               continue;
           end
        end
        
        if ExcludeIfHasSameDirBeforeDiffDir==1
           if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==5;
               % then throw out
                 disp(['SKIPPING ' birdname '- ' exptname ': DID SAME DIR BEFORE BIDIR']);
                continue;
            end
        end
        
        if ExcludeNotFullyLabeled==1
            if strcmp(birdname, 'bk34bk68') && strcmp(exptname, 'SeqDepPitchLMAN3');
                disp(['SKIPPING ' birdname '- ' exptname ': NOT FULLY LABELED YET learning']); % reprobed during trajectory, day 3.
                continue;
            end
            if strcmp(birdname, 'gr41gr90') && strcmp(exptname, 'SeqDepPitchLMAN2');
                disp(['SKIPPING ' birdname '- ' exptname ': NOT FULLY LABELED YET learning']); % catch left high at 0.8 on day 4, could potentailly use for day 5. 
                continue;
            end
            
        end
        
        if ExcludeIfFirstTargDriveMore==1
             if strcmp(birdname, 'pu64bk13') && strcmp(exptname, 'SeqDepPitchShift');
                disp(['SKIPPING ' birdname '- ' exptname ': FIRST TARG DROVE MORE...']); % reprobed during trajectory, day 3.
                continue;
            end
        end
        
            
        
        % What days to plot? (start, N days before start multidir; end:
        % N days at end of multidir)
        DayBeforeMDon_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        MDoff_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
        
        daystokeep=DayBeforeMDon_ind-DaysToPlot(1)+1:DayBeforeMDon_ind+DaysToPlot(2); % will take N days BEFORE start. end at last day (not post).
        
        % trim end if reaches end of bidir
        if daystokeep(end)>MDoff_ind;
            if UseOldVersion==1
                daystokeep(end) = MDoff_ind;
            else
                daystokeep= daystokeep(daystokeep<=MDoff_ind);
            end
        end
        
        if ThrowOutIfNotEnoughData==1
        if length(daystokeep)<sum(DaysToPlot)
            continue
        end
        end
        
        % === save the number of WN single dir days that occurred before
        % bidir on
        WNDay1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        NumSingleTargDays=DayBeforeMDon_ind-WNDay1+1;
        
        % -- collect
        
        
        % ================== COLLECT
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        multidirsyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targLearnDir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        if strcmp(multidir_syls{1}, targsyl)
            secondtarg=multidir_syls{2};
        elseif strcmp(multidir_syls{2}, targsyl)
            secondtarg=multidir_syls{1};
        else
            disp('PROBLEM - which targ is secondtarg?');
            dasfacas
        end
        
        % --- if is LMAN, use within time field
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            FFdataField='meanFF_DevFromBase_WithinTimeWindow';
        else
            FFdataField='meanFF_DevFromBase';
        end
        
        
        % == TARG1
        if max(daystokeep)>length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).(FFdataField));
            daystokeep=daystokeep(daystokeep<=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).(FFdataField)));
            disp([birdname '-' exptname]);
        end
        
%         if strcmp(birdname, 'pu53wh88')
%             keyboard
%         end
        FFvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).(FFdataField)(daystokeep);
        
        
        % --- make sure is one column (LMAN and not LMAN differ)
        tmp=size(FFvals);
        if tmp(1)~=tmp(2)
            if tmp(2)>1;
                FFvals=FFvals';
            end
        end
        FFvals=FFvals.*targLearnDir;
        
        DATSTRUCT.data.FirstTarg(Exptcount).MeanFFRelBase=FFvals;
        
        % VALUE TO NORM ALL TO
        ValueToNormAllTo=FFvals(DaysToPlot(1));
        FFvals=FFvals./ValueToNormAllTo;
        DATSTRUCT.data.FirstTarg(Exptcount).MeanFFRelBase_NormTargLastBaseDay=FFvals;
        
        % ------ COUNT SAMPLE SIZE (SEQUENCE LEARNING?)
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            Numrends_syl1_end=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{daystokeep(end)});
            Numrends_syl1_start=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow{daystokeep(1)});
        else
            Numrends_syl1_end=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals{daystokeep(end)});
            Numrends_syl1_start=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals{daystokeep(1)});
        end
        
        
        % == TARG2
        FFvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).(FFdataField)(daystokeep);
        % --- make sure is one column (LMAN and not LMAN differ)
        tmp=size(FFvals);
        if tmp(1)~=tmp(2)
            if tmp(2)>1;
                FFvals=FFvals';
            end
        end
        FFvals=FFvals.*targLearnDir;

        DATSTRUCT.data.SecondTarg(Exptcount).MeanFFRelBase=FFvals;
        % signifincatly diff from 0?
        ffvalsRaw = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).FFvals(daystokeep);
        ffbase = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(secondtarg).meanFF_RecalcDays;

        pvalsAll = [];
        FFsemAll = [];
        for k=1:length(ffvalsRaw)
            ffvalstmp = cell2mat(ffvalsRaw{k});
            ffvalstmp = ffvalstmp - ffbase;
        if isempty(ffvalstmp)
            p = nan;
            FFsemAll = [FFsemAll nan];
            continue
        end
    
            p = signrank(ffvalstmp);
            pvalsAll = [pvalsAll p];
            FFsemAll = [FFsemAll lt_sem(ffvalstmp)];
        end
        DATSTRUCT.data.SecondTarg(Exptcount).pvalsAll_signrankRelZero=pvalsAll;
        DATSTRUCT.data.SecondTarg(Exptcount).FFsemAll=FFsemAll;
            
        
        
        % VALUE TO NORM ALL TO
        FFvals=FFvals./ValueToNormAllTo;
        DATSTRUCT.data.SecondTarg(Exptcount).MeanFFRelBase_NormTargLastBaseDay=FFvals;
        
         % ------ COUNT SAMPLE SIZE (SEQUENCE LEARNING?)
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            Numrends_syl2_end=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).FFvals_WithinTimeWindow{daystokeep(end)});
            Numrends_syl2_start=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).FFvals_WithinTimeWindow{daystokeep(1)});
        else
            Numrends_syl2_end=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).FFvals{daystokeep(end)});
            Numrends_syl2_start=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(secondtarg).FFvals{daystokeep(1)});
        end
       
        disp(' ');
        disp([birdname '-' exptname]);
        disp(['START (syl1/syl2): ' num2str(Numrends_syl1_start) '/' num2str(Numrends_syl2_start) ' --> ' num2str(Numrends_syl1_start/(Numrends_syl1_start+Numrends_syl2_start))]);
        disp(['END (syl1/syl2)  : ' num2str(Numrends_syl1_end) '/' num2str(Numrends_syl2_end) ' --> ' num2str(Numrends_syl1_end/(Numrends_syl1_end+Numrends_syl2_end))]);

        
        
        % == OTHERS (SAME)
        FFvals_All=[];
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            % skip if is targ1 or 2
            if strcmp(syl, targsyl) | strcmp(syl, secondtarg);
                continue
            end
            
            % confirm it is same type
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                continue
            end
            %                         if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab==0;
            %                 continue
            %             end
            
            % ----- Collect data
            FFvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(FFdataField)(daystokeep);
            % --- make sure is one column (LMAN and not LMAN differ)
            tmp=size(FFvals);
            if tmp(1)~=tmp(2)
                if tmp(2)>1;
                    FFvals=FFvals';
                end
            end
            
            FFvals_All=[FFvals_All FFvals];
        end
        FFvals_All=FFvals_All.*targLearnDir;
        DATSTRUCT.data.OthersSameType(Exptcount).MeanFFRelBase=nanmean(FFvals_All, 2);
        
        if adHocPlotSameType==1

            SameType_Others = [SameType_Others FFvals_All];
        ExptcountAll = [ExptcountAll Exptcount*ones(1,size(FFvals_All,2))];
        end
        
        % VALUE TO NORM ALL TO
        FFvals=nanmean(FFvals_All, 2)./ValueToNormAllTo;
        DATSTRUCT.data.OthersSameType(Exptcount).MeanFFRelBase_NormTargLastBaseDay=FFvals;
        
        
        % == OTHERS (DIFF)
        FFvals_All=[];
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            % skip if is targ1 or 2
            if strcmp(syl, targsyl) | strcmp(syl, secondtarg);
                continue
            end
            
            % confirm it is diff type
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1;
                continue
            end
            
            %             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab==1;
            %                 continue
            %             end
            % ----- Collect data
            FFvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(FFdataField)(daystokeep);
            % --- make sure is one column (LMAN and not LMAN differ)
            tmp=size(FFvals);
            if tmp(1)~=tmp(2)
                if tmp(2)>1;
                    FFvals=FFvals';
                end
            end
            
            FFvals_All=[FFvals_All FFvals];
        end
        FFvals_All=FFvals_All.*targLearnDir;
        DATSTRUCT.data.OthersDiffType(Exptcount).MeanFFRelBase=nanmean(FFvals_All, 2);
        
                if adHocPlotSameType==0
                       
        SameType_Others = [SameType_Others FFvals_All];
        ExptcountAll = [ExptcountAll Exptcount*ones(1,size(FFvals_All,2))];
                end

        
        % VALUE TO NORM ALL TO
        FFvals=nanmean(FFvals_All, 2)./ValueToNormAllTo;
        DATSTRUCT.data.OthersDiffType(Exptcount).MeanFFRelBase_NormTargLastBaseDay=FFvals;
        
        
        
        
        % == OTHERS (ALL)
        FFvals_All=[];
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            % skip if is targ1 or 2
            if strcmp(syl, targsyl) | strcmp(syl, secondtarg);
                continue
            end
            
            % ----- Collect data
            FFvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(FFdataField)(daystokeep);
            % --- make sure is one column (LMAN and not LMAN differ)
            tmp=size(FFvals);
            if tmp(1)~=tmp(2)
                if tmp(2)>1;
                    FFvals=FFvals';
                end
            end
            
            FFvals_All=[FFvals_All FFvals];
        end
        FFvals_All=FFvals_All.*targLearnDir;
        DATSTRUCT.data.OthersAll(Exptcount).MeanFFRelBase=nanmean(FFvals_All, 2);
        
        % VALUE TO NORM ALL TO
        FFvals=nanmean(FFvals_All, 2)./ValueToNormAllTo;
        DATSTRUCT.data.OthersAll(Exptcount).MeanFFRelBase_NormTargLastBaseDay=FFvals;
        
        
        
        % ===== INFORMATION
        DATSTRUCT.information(Exptcount).birdname=birdname;
        DATSTRUCT.information(Exptcount).exptname=exptname;
        DATSTRUCT.information(Exptcount).dayIndsUsed=daystokeep;
        DATSTRUCT.information(Exptcount).numPreBidirDays=DaysToPlot(1);
        DATSTRUCT.information(Exptcount).numBidirDays=DaysToPlot(2);
        DATSTRUCT.information(Exptcount).targLearnDir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        DATSTRUCT.information(Exptcount).NumSingleTargDaysBeforeBidir=NumSingleTargDays;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(secondtarg).similar_to_targ==1;
        DATSTRUCT.information(Exptcount).targ2_sametype_rel_targ1=1;
        else
        DATSTRUCT.information(Exptcount).targ2_sametype_rel_targ1=0;
        end
        
        presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(secondtarg).presyl_similar_to_targ_presyl;
        DATSTRUCT.information(Exptcount).targ2_presim_rel_targ1=presim;
        
        
        durThrAdaptingSingleCtxt=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumDaysAdapticeThr_col1OK_col2complete(2);
        DATSTRUCT.information(Exptcount).durThrAdaptingSingleCtxt=durThrAdaptingSingleCtxt;
        
        
        
        % =====
        
        Exptcount=Exptcount+1;
    end
end


%% ==== PLOT DATAPOINTS FOR EACH DAY (BARS, REL TO BIDIR START)

numexpts=length(DATSTRUCT.data.FirstTarg);
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


for i=1:numexpts
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([DATSTRUCT.information(i).birdname '-' DATSTRUCT.information(i).exptname]);
    
    X=DATSTRUCT.information(i).dayIndsUsed;
    
    if ~DATSTRUCT.information(i).targ2_sametype_rel_targ1
        lt_plot_text(0, 0.5, 'skipped, since not sametype', 'r');
        continue
    end
    
    % ----- targ 1
    Y=DATSTRUCT.data.FirstTarg(i).MeanFFRelBase;
    
    lt_plot(X, Y, {'Color', 'k','LineStyle','-'});
    
    % ---- targ 2
    Y=DATSTRUCT.data.SecondTarg(i).MeanFFRelBase;
    
    lt_plot(X, Y, {'Color', 'k','LineStyle','--'});
    
    % --- others (same);
    Y=DATSTRUCT.data.OthersSameType(i).MeanFFRelBase;
    if ~isempty(Y);
        plot(X, Y, 'o-b');
    end
    
    % --- others (siff);
    Y=DATSTRUCT.data.OthersDiffType(i).MeanFFRelBase;
    if ~isempty(Y);
        plot(X, Y, 'o-r');
    end
    
    % ---
    line(xlim, [0 0], 'Color','m');
    
    % ---- line for bidir start
    xx=DATSTRUCT.information(i).numPreBidirDays+DATSTRUCT.information(i).dayIndsUsed(1)-1;
    
    line([xx+0.5 xx+0.5], ylim);
    
end


%% ==== PLOT DATAPOINTS FOR EACH DAY (BARS, REL TO BIDIR START) [OVERLAY ALL EXPERIMENTS
numexpts=length(DATSTRUCT.data.FirstTarg);

lt_figure; hold on;
onlyKeepExptIfFullData =1;
if onlyKeepExptIfFullData==1
    title('only expt with at least this much data');
else
    title('all expt, regardless if enough data');
end
xlabel('day (bidir start at 1)');
ylabel('change in FF, rel baseline (red means p<0.005 from base)');


Yall_targ = [];
Yall_sametype = [];

for i=1:numexpts
    X=DATSTRUCT.information(i).dayIndsUsed;
    bidirStart = DATSTRUCT.information(i).numPreBidirDays+DATSTRUCT.information(i).dayIndsUsed(1)-1;
    X = X - bidirStart;
    
    if ~DATSTRUCT.information(i).targ2_sametype_rel_targ1
        continue
    end
    
    if onlyKeepExptIfFullData==1
        
        if sum(~isnan(DATSTRUCT.data.FirstTarg(i).MeanFFRelBase))<sum(DaysToPlot)
            continue
        end
    disp(DATSTRUCT.data.FirstTarg(i).MeanFFRelBase)
    end
    
    
    % ----- targ 1
    Y=DATSTRUCT.data.FirstTarg(i).MeanFFRelBase;
    
    plot(X, Y, '-k', 'LineWidth', 1.5);
    Yall_targ = [Yall_targ; Y'];
    
    
    
    % ---- targ 2
    Y=DATSTRUCT.data.SecondTarg(i).MeanFFRelBase;
    
    pValsAll = DATSTRUCT.data.SecondTarg(i).pvalsAll_signrankRelZero;
    FFsemAll = DATSTRUCT.data.SecondTarg(i).FFsemAll;
    
    indstmp = pValsAll<0.005; % those that are diff from baseline (0))\
    if sum(indstmp) >0;
        lt_plot(X(indstmp), Y(indstmp), {'Errors', FFsemAll(indstmp), 'Color','r'});
    end
    
%     indstmp = pValsAll>=0.05;
%     if sum(indstmp) >0;
%         lt_plot(X(indstmp), Y(indstmp), {'Errors', FFsemAll(indstmp), 'MarkerFaceColor', 'none'});
%     end
    
    
    plot(X, Y, '-b', 'LineWidth', 2);
    Yall_sametype = [Yall_sametype; Y'];
    disp(DATSTRUCT.information(i).dayIndsUsed)
    %     % --- others (same);
    %     Y=DATSTRUCT.data.OthersSameType(i).MeanFFRelBase;
    %     if ~isempty(Y);
    %         plot(X, Y, '-m');
    %     end
    %
    %     % --- others (siff);
    %     Y=DATSTRUCT.data.OthersDiffType(i).MeanFFRelBase;
    %     if ~isempty(Y);
    %         plot(X, Y, '-r');
    %     end
end

% === plot mean
% lt_plot(X, mean(Yall_targ), {'Errors', lt_sem(Yall_targ)});
% lt_plot(X, mean(Yall_sametype), {'Errors', lt_sem(Yall_sametype)});

% ---
axis normal
line(xlim, [0 0], 'Color','m');
line([0.5 0.5], ylim, 'Color', 'm');


% ======================= PLOT AS BARS ACROSS DAYS
lt_figure; hold on;

% - targ
lt_subplot(2,1,1); hold on;
lt_plot_bar(1:size(Yall_targ,2), mean(Yall_targ), {'Errors', lt_sem(Yall_targ)})
ylim([-200 250]);

% - nontarg
lt_subplot(2,1,2); hold on;
lt_plot_bar(1:size(Yall_sametype,2), mean(Yall_sametype), {'Errors', lt_sem(Yall_sametype)})
ylim([-200 250]);
% -- which diff from 0?
[~, p] = ttest(Yall_sametype)
% -- compare 


% ======================= COMPARE 2 DAYS
if DaysToPlot(2)>=7
daytmp1 = -1;
daytmp2 = 10;

% -- add num base days
daytmp1 = daytmp1 + DaysToPlot(1);
daytmp2 = daytmp2 + DaysToPlot(1);

lt_figure; hold on;
YdiffAllTmp = {}; % to collect differences

% --- targ
X = [1 2];

ytmp = Yall_targ(:, [daytmp1 daytmp2]);
[h, p] = ttest(ytmp(:,1), ytmp(:,2))
plot(X, ytmp, '-ok');
lt_plot_bar(X, mean(ytmp), {'Errors', lt_sem(ytmp)});
lt_plot_pvalue(p, 'ttest(targ)', 1);

YdiffAllTmp{1} = diff(ytmp, 1, 2);


% --- same type
X = [4 5];
ytmp = Yall_sametype(:, [daytmp1 daytmp2]);
[h, p] = ttest(ytmp(:,1), ytmp(:,2))
plot(X, ytmp, '-ok');
lt_plot_bar(X, mean(ytmp), {'Errors', lt_sem(ytmp)});
lt_plot_pvalue(p, 'ttest(same)', 1);
% -- which same types diff from 0?


else
    disp('NOT PLOTTING FURTHER LEARNING!! (supplemental fig)');
end
%% ==== PLOT DATAPOINTS FOR EACH DAY (BARS, REL TO BIDIR START) [USING NORM TO TARG]

numexpts=length(DATSTRUCT.data.FirstTarg);
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];



for i=1:numexpts
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([DATSTRUCT.information(i).birdname '-' DATSTRUCT.information(i).exptname]);
    
    X=DATSTRUCT.information(i).dayIndsUsed;
    
    % ----- targ 1
    Y=DATSTRUCT.data.FirstTarg(i).MeanFFRelBase_NormTargLastBaseDay;
    
    lt_plot(X, Y, {'Color', 'k','LineStyle','-'});
    
    % ---- targ 2
    Y=DATSTRUCT.data.SecondTarg(i).MeanFFRelBase_NormTargLastBaseDay;
    
    lt_plot(X, Y, {'Color', 'k','LineStyle','--'});
    
    % --- others (same);
    Y=DATSTRUCT.data.OthersSameType(i).MeanFFRelBase_NormTargLastBaseDay;
    if ~isempty(Y);
        plot(X, Y, 'o-b');
    end
    
    % --- others (siff);
    Y=DATSTRUCT.data.OthersDiffType(i).MeanFFRelBase_NormTargLastBaseDay;
    if ~isempty(Y);
        plot(X, Y, 'o-r');
    end
    
    % ---
    line(xlim, [0 0], 'Color','m');
    
    % ---- line for bidir start
    xx=DATSTRUCT.information(i).numPreBidirDays+DATSTRUCT.information(i).dayIndsUsed(1)-1;
    
    line([xx+0.5 xx+0.5], ylim);
    
    
end


%% ==== DISTRIBUTION OF SINGLE TARG PHASE NUM DAYS
lt_figure; hold on;

inds=find([DATSTRUCT.information.targ2_sametype_rel_targ1]);

% -- num days of single context before start dual
lt_subplot(3,2,1);
hold on;
title('num single targ days before start bidir [same type pairs only');
Y=[DATSTRUCT.information(inds).NumSingleTargDaysBeforeBidir];

lt_plot_histogram(Y);
lt_plot_annotation(1, ['mean=' num2str(mean(Y)) '; std=' num2str(std(Y)) '; min= ' num2str(min(Y)) '; max=' num2str(max(Y))], 'k')

lt_plot_annotation(4, ['n=' num2str(length(inds))]);

% ======== num days single context was adaptive changing of thr
lt_subplot(3,2,2);
hold on;
title('num single targ days adaptive change of thr');
Y=[DATSTRUCT.information(inds).durThrAdaptingSingleCtxt];
lt_plot_histogram(Y);
lt_plot_annotation(1, ['mean=' num2str(mean(Y)) '; std=' num2str(std(Y)) '; min= ' num2str(min(Y)) '; max=' num2str(max(Y))], 'k')


% --- display names of all expts used for dual context (i.e. same type)
inds=find([DATSTRUCT.information.targ2_sametype_rel_targ1]);
for ind=inds
    
    bname=DATSTRUCT.information(ind).birdname;
    ename=DATSTRUCT.information(ind).exptname;
    
    disp([bname '-' ename])
end


%% +++++++++++++ USING RAW FF
%% ====== ALL

NormToTarg=0;
PlotSameTypeOnly=0;


lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1
lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'All - raw FF', 'k')




%% ====== SAME TYPE

NormToTarg=0;
PlotSameTypeOnly=1; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1
lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'Same types - raw FF', 'k')

%% ====== SAME TYPE, trying different ways to plot

NormToTarg=0;
PlotSameTypeOnly=1; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1_1
lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'Same types - raw FF', 'k')


%% ==== DIFF TYPE

NormToTarg=0;
PlotSameTypeOnly=2; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1

lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'Diff types - raw FF', 'k')

%% ====== SAME TYPE SAME SEQ

NormToTarg=0;
PlotSameTypeOnly=3; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1

lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'SameType SameSeq - raw FF', 'k')


%% ====== SAME TYPE DIFF SEQ

NormToTarg=0;
PlotSameTypeOnly=4; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1

lt_subplot(3,3,9); hold on;
lt_plot_annotation(1, 'SameType DiffSeq - raw FF', 'k')



%% +++++++++++++++++++++++++++ USING NORM FF
%% ====== COMBINE ALL EXPERIMENTS (LOCKED TO BIDIR START) [norm to targ 1 last bline day]

%% ====== ALL

NormToTarg=1;
PlotSameTypeOnly=0;


lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1
lt_subtitle('All - norm FF');




%% ====== SAME TYPE

NormToTarg=1;
PlotSameTypeOnly=1; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1
lt_subtitle('Same types - norm FF');

%% ==== DIFF TYPE

NormToTarg=1;
PlotSameTypeOnly=2; % 0 = all; 1= same type; 2= diff type

lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR_v2sub1
lt_subtitle('Diff types - norm FF');



%% === PLOT SEPARATION 
% -- separation on last single targ day vs. separation on last bidir day
% === PARAMS
PlotSameTypeOnly = 1;
lt_figure; hold on;
title('separation, last day of single and dual');

% === RUNS
if PlotSameTypeOnly==1;
    ExptInds=find([DATSTRUCT.information.targ2_sametype_rel_targ1]);
else
    disp('PROBLEM, DO WHAT?');
end


lastSingleDay = DATSTRUCT.information(1).numPreBidirDays;
lastBidirDay = DATSTRUCT.information(1).numBidirDays + DATSTRUCT.information(1).numPreBidirDays;

Y = nan(length(ExptInds), 2);

% ==== end of single targ phase
dayind = lastSingleDay;
X = 1;

dattmp1 = [DATSTRUCT.data.FirstTarg(ExptInds).MeanFFRelBase];
dattmp1 = dattmp1(dayind, :)';

dattmp2 = [DATSTRUCT.data.SecondTarg(ExptInds).MeanFFRelBase];
dattmp2 = dattmp2(dayind, :)';

Y(:, X) = dattmp1 - dattmp2;

% ==== end of dual targ phase
dayind = lastBidirDay;
X = 2;

dattmp1 = [DATSTRUCT.data.FirstTarg(ExptInds).MeanFFRelBase];
dattmp1 = dattmp1(dayind, :)';

dattmp2 = [DATSTRUCT.data.SecondTarg(ExptInds).MeanFFRelBase];
dattmp2 = dattmp2(dayind, :)';

Y(:, X) = dattmp1 - dattmp2;

% === plot
plot([1 2], Y, '-k')

lt_plot_bar([1 2], mean(Y,1), {'Errors', lt_sem(Y)})

% === stats
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'signrank', 1);
lt_plot_text(1, 1.1*max(max(Y)), [num2str(mean(Y(:,1))) '(' num2str(lt_sem(Y(:,1))) ')'], 'b')
lt_plot_text(2, 1.1*max(max(Y)), [num2str(mean(Y(:,2))) '(' num2str(lt_sem(Y(:,2))) ')'], 'b')


% BELOW: actuall comparing last day of isngle and dual separately for each
% syl type ...
% % ==== targ
% X = [1 2];
% plotcol = 'k';
% targfield = 'FirstTarg';
% 
% % --- 
% datarray = [DATSTRUCT.data.(targfield)(ExptInds).MeanFFRelBase]; % days x expt
% datarray = datarray([lastSingleDay lastBidirDay], :)';
% 
% plot(X, datarray, '-o', 'Color', plotcol);
% 
% % ==== same
% X = [4 5];
% plotcol = 'b';
% targfield = 'SecondTarg';
% 
% % --- 
% datarray = [DATSTRUCT.data.(targfield)(ExptInds).MeanFFRelBase]; % days x expt
% datarray = datarray([lastSingleDay lastBidirDay], :)';
% 
% plot(X, datarray, '-o', 'Color', plotcol);
% 







%% ===== CORRELATION BETWEEN TARG 2 SHIFT AND TARG 1 SHIFT
% === for each day, plot shift from yesterday, correalted between syls?


figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

PlotSameTypeOnly=0;
if PlotSameTypeOnly==1;
tmp=[DATSTRUCT.information.targ2_sametype_rel_targ1];
ExptInds=find(tmp==1);
elseif PlotSameTypeOnly==2;
    tmp=[DATSTRUCT.information.targ2_sametype_rel_targ1];
ExptInds=find(tmp==0);
else
    ExptInds=1:length(DATSTRUCT.information);
end
    

Y_DayDiffs_FirstTarg=diff([DATSTRUCT.data.FirstTarg(ExptInds).MeanFFRelBase]); % [pitch diffs x expt];
Y_DayDiffs_SecondTarg=diff([DATSTRUCT.data.SecondTarg(ExptInds).MeanFFRelBase]);


for day = 1:length(DATSTRUCT.data.FirstTarg(1).MeanFFRelBase);
    if day==1;
        continue
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['day ' num2str(day) ' minus ' num2str(day-1)]);
    xlabel('targ 2 shift');
    ylabel('targ 1 shift');
    
    
    % shift from yesterday
    Y_firsttarg=Y_DayDiffs_FirstTarg(day-1,:);
    Y_secondtarg=Y_DayDiffs_SecondTarg(day-1, :);
    
    lt_regress(Y_firsttarg, Y_secondtarg, 1, 0, 1, 1, 'k');
    lt_plot_45degScatter(Y_secondtarg, Y_firsttarg);    
end

    lt_subtitle('day to day shifts enhanced during WN?');
    
    
 
    %% == adhoc for revisions [SAME TYPE, NON TARG] 
    % breaking out into individual syl (otherwise exactly same as
    % same-type, nontarg in DATSTRUCT. In latter each expt taken mean.
    
    tmp=[DATSTRUCT.information.targ2_sametype_rel_targ1];
ExptInds=find(tmp==1);

inds = ismember(ExptcountAll, ExptInds);
Y = SameType_Others(:,inds);
lt_figure; hold on; 
if adHocPlotSameType==1
    title('same type, nontarg, each syl')
else
        title('diff type, nontarg, each syl')
end

lt_subplot(2,1,1); hold on;
plot(1:size(Y,1), Y, '-k');
Ymean = mean(Y,2);
lt_plot_bar(1:length(Ymean), Ymean', {'Errors', lt_sem(Y')})

lt_subplot(2,1,2); hold on;
Ymean = mean(Y,2);
lt_plot_bar(1:length(Ymean), Ymean', {'Errors', lt_sem(Y')})





