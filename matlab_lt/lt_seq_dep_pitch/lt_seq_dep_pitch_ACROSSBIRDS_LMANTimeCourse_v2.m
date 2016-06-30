function lt_seq_dep_pitch_ACROSSBIRDS_LMANTimeCourse_v2(SeqDepPitch_AcrossBirds, Params, DayToLockAllTo)
%% LT 12/22/15

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% DISPLAY STATS FOR WHEN THERE WERE INACTIVATION DAYS FOR ALL BIRDS
disp(' ');
disp('WN day 1  --  Days with musc (WNday1=1)  --  Bidir day 1(WNday1=1)');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        bidir_day1=nan;
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
            bidir_day1=bidir_day1-WNday1+1;
        end
        
        
        disp([birdname '-' exptname ':  ' num2str(WNday1) ' -- ' mat2str(days_with_musc_relWNday1) ' -- ' num2str(bidir_day1)]);
        
    end
end


%% PLOT TIMELINE FOR ALL EXPERIMENTS [IN PROGRESS, ADDING REVERSION AT TARGET, SEE LINE WITH: FFminusBase_MUSC]
% EACH ROW IS ONE EXPERIMENT, ALL LOCKED TO WN ON. SHOW EACH INACTIVATION
% DAY. SHOW BIDIR. SHOW CV DROP AND REVERSION AT TARGET, AND PERCENT MAX
% LEARNING BY TARGET

% WN day 1 =  ind1

lt_figure; hold on;
BirdnamesExptname_all={};
title('Timeline (muscimol in bold) (bidir in blue)');
xlabel('day (WN day 1 = 1)');



CountExpt=1;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        days_with_musc_relWNday1=days_with_musc-WNday1+1;
        
        % --- bidir start and end
        bidir_day1=nan;
        bidir_lastday=nan;
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
            bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
            bidir_day1=bidir_day1-WNday1+1;
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_end_Inds');
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds(end);
            else
                bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            bidir_lastday=bidir_lastday-WNday1+1;
        end
        
        % -- last day
        if ~isnan(bidir_lastday);
            WN_last=bidir_lastday;
        else
            WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            WN_last=WN_last-WNday1+1;
        end
        
        % ---- first day (relative to WN on day)
        firstday=-(WNday1-2);
        
        
        % ======= REVERSION AT TARGET
        FFminusBase_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        FFminusBase_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        
        
        % === plot
        % - line for all days with any data
        X=firstday:WN_last;
        Y=CountExpt*ones(1,length(X));
        
        plot(X, Y, '-ok', 'LineWidth', 1);
        
        
        
        % -- line for start and end of bidir
        line([bidir_day1-0.5 bidir_day1-0.5], [CountExpt-0.4 CountExpt+0.4]);
        line([bidir_lastday+0.5 bidir_lastday+0.5], [CountExpt-0.4 CountExpt+0.4]);
        
        
        % -- fill in musc days
        X=days_with_musc_relWNday1;
        Y=CountExpt*ones(1, 1:length(X));
        
        plot(X, Y, 'ok', 'MarkerFaceColor','k')
        
        % --- line for end of WN
        line([WN_last+0.5 WN_last+0.5], [CountExpt-0.4 CountExpt+0.4], 'Color','r');
        
        
        % -- collect bird name and expt
        BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
        
        CountExpt=CountExpt+1;
    end
end

set(gca, 'YTick', 1:CountExpt-1);
set(gca, 'YTickLabel', BirdnamesExptname_all);


% ---- WN on line
line([0.5 0.5], ylim, 'Color','r');



%%  for all days plot shift + reversion at target [THINKING OF PUTTING IN CODE ABOVE INSTEAD. SEE ABOVE, STILL IN PROGRESS]
if (0)
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        
        for ii=1:numexpts;
            
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            days_with_musc_relWNday1=days_with_musc-WNday1+1;
            
            % --- bidir start and end
            bidir_day1=nan;
            bidir_lastday=nan;
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
                bidir_day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
                bidir_day1=bidir_day1-WNday1+1;
                
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_end_Inds');
                    bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds(end);
                else
                    bidir_lastday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
                end
                bidir_lastday=bidir_lastday-WNday1+1;
            end
            
            % -- last day
            if ~isnan(bidir_lastday);
                WN_last=bidir_lastday;
            else
                WN_last=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
                WN_last=WN_last-WNday1+1;
            end
            
            % ---- first day (relative to WN on day)
            firstday=-(WNday1-2);
            
            % -- collect bird name and expt
            BirdnamesExptname_all=[BirdnamesExptname_all [birdname(1:4) '-' exptname(end-4:end)]];
            
            CountExpt=CountExpt+1;
        end
    end
    
end


%% ===== GET ALL DATA LOCKED TO SPECIFIC DAY

disp(' ----- ')
% disp(' BIDIR CONSOLIDATION ');
%
ExptCount=1;


DATSTRUCT=struct;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        % === GET DAY INDS RELATIVE TO A LOCK DAY (DAY1)
        
        if strcmp(DayToLockAllTo, 'BidirDay1')
            
            % 1) skip bird if no bidir or if more than two targ
            if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
                disp([birdname '-' exptname '- SKIPPED bidir consolidation (no bidir data)']);
                continue
            end
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.num_multidir_syls~=2;
                disp([birdname '-' exptname '- SKIPPED num multidir syls is not 2']);
                continue
            end
            
            % 2) get dates
            FirstDay_AllElseLockedToThis=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1);
            
        elseif strcmp(DayToLockAllTo, 'WNday1');
            
            FirstDay_AllElseLockedToThis=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            
        elseif strcmp(DayToLockAllTo, 'SingleDirConsolid');
            
            % --- need to enter consolid dates in params
            birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
            exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
            if ~isempty(ind);
                consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
                
            else
                % skip bird
                disp(['skipped ' birdname '-' exptname ' - NO CONSOLID PERIOD ENTERED']);
                continue
            end
            
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            FirstDay_AllElseLockedToThis=WNday1+consolPeriod(1)-1;
        
        elseif strcmp(DayToLockAllTo, 'BaseDay1');
            
            FirstDay_AllElseLockedToThis=1;
        end
        
        
        
        % ---- GET OTHER DATES
        LastDayWN_ActualInds=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow);
        LastDayWN_RelInds=LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1;
        
        
        
        % ===== COLLECT DATA [FIRST TARGET
        
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean_PBS=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean_MUSC=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        
        % ---- PBS
        ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);
        
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean_PBS=ffmeans_RelInds;
        
        % ---- MUSC
        ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow;
        ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);
        
        DATSTRUCT.firsttarget(ExptCount).FFminusBase_Mean_MUSC=ffmeans_RelInds;
        
        

        % ====== NONTARGETS (PUT INTO CLASSES) ============================
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        SAMETYPEpbs=[];
        SAMETYPEmusc=[];
        DIFFTYPEpbs=[];
        DIFFTYPEmusc=[];
        ALLpbs=[];
        ALLmusc=[];
        
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean_PBS=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean_PBS=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean_PBS=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean_MUSC=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean_MUSC=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean_MUSC=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);

       
        for k=1:length(SylsUnique);
            syl=SylsUnique{k};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                continue;
            end
            
            % ====== [PBS] extract ff information
            ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow;
            ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);
            
            
            % --- SAME TYPE
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                SAMETYPEpbs=[SAMETYPEpbs; ffmeans_RelInds];
                % --- DIFF TYPE
            else SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                DIFFTYPEpbs=[DIFFTYPEpbs; ffmeans_RelInds];
            end
            % --- ALL
            ALLpbs=[ALLpbs; ffmeans_RelInds];
            
            
            % ====== [MUSC] extract ff information
            ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow;
            ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);

            % --- SAME TYPE
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                SAMETYPEmusc=[SAMETYPEmusc; ffmeans_RelInds];
                % --- DIFF TYPE
            else SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                DIFFTYPEmusc=[DIFFTYPEmusc; ffmeans_RelInds];
            end
            % --- ALL
            ALLmusc=[ALLmusc; ffmeans_RelInds];
        end
        
        % ==== OUTPUT - nontargets =================================
        % get mean and sem
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean_PBS=nanmean(SAMETYPEpbs, 1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean_PBS=nanmean(DIFFTYPEpbs, 1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean_PBS=nanmean(ALLpbs, 1);
        
        DATSTRUCT.nontargets_Same(ExptCount).FFminusBase_Mean_MUSC=nanmean(SAMETYPEmusc, 1);
        DATSTRUCT.nontargets_Diff(ExptCount).FFminusBase_Mean_MUSC=nanmean(DIFFTYPEmusc, 1);
        DATSTRUCT.nontargets_All(ExptCount).FFminusBase_Mean_MUSC=nanmean(ALLmusc,1 );
        
        
        
        
        
        
        % ==== IF WANT BIDIR, THEN GET SECOND TARGET TOO
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds') ...
                & SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.num_multidir_syls==2;
            
            
            %  --- what is other target?
            syl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls{1};
            syl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls{2};
            
            if strcmp(syl1, targsyl)
                othersyl=syl2;
            elseif strcmp(syl2,targsyl)
                othersyl=syl1;
            else
                disp('PROBLEM!!!');
                dafascawr3a;
            end
            
            
            % ---------------- COLLECT
            DATSTRUCT.secondtarg(ExptCount).FFminusBase_Mean_PBS=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
            DATSTRUCT.secondtarg(ExptCount).FFminusBase_Mean_MUSC=nan(1, LastDayWN_ActualInds-FirstDay_AllElseLockedToThis+1);
            
            % ---- PBS
            ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(othersyl).meanFF_DevFromBase_WithinTimeWindow;
            ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);
            
            DATSTRUCT.secondtarget(ExptCount).FFminusBase_Mean_PBS=ffmeans_RelInds;
            
            % ---- MUSC
            ffmeans_actualInds=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(othersyl).meanFF_DevFromBase_WithinTimeWindow;
            ffmeans_RelInds=ffmeans_actualInds(FirstDay_AllElseLockedToThis:LastDayWN_ActualInds);
            
            DATSTRUCT.secondtarget(ExptCount).FFminusBase_Mean_MUSC=ffmeans_RelInds;
            
        end
        
        
       
        
        
        
        % ====== SAVE INFORMATION ABOUT THIS EXPERIMENT
        
        DATSTRUCT.INFORMATION(ExptCount).birdname=birdname;
        DATSTRUCT.INFORMATION(ExptCount).exptname=exptname;
        
        % BIDIR END, OR WN END (IF NO BIDIR)
        LastDayWN_ActualInds_actual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        if  strcmp(DayToLockAllTo, 'BidirDay1')
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_end_Inds');
                LastDayWN_ActualInds_actual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds;
                %         disp([ num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds) ' - ' num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd)]);
                
            end
        end
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
            BidirStartActual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds;
            %         disp([ num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds) ' - ' num2str(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd)]);
            BidirStart_Rel=BidirStartActual-FirstDay_AllElseLockedToThis+1;
            DATSTRUCT.INFORMATION(ExptCount).BidirDay1_RelInd=BidirStart_Rel;
        end
        
        
        
        DATSTRUCT.INFORMATION(ExptCount).LastWNDay_RelInds=LastDayWN_RelInds;
        DATSTRUCT.INFORMATION(ExptCount).LastDay_RelInds_BidirIfExist_OtherwiseWN=LastDayWN_ActualInds_actual-FirstDay_AllElseLockedToThis+1;
        DATSTRUCT.INFORMATION(ExptCount).targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        
        ExptCount=ExptCount+1;
        
    end
end
%% ====== PLOT [targ and second targ]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT.firsttarget);
    
    birdname=DATSTRUCT.INFORMATION(i).birdname;
    exptname=DATSTRUCT.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day');
    ylabel('shift from baseline (hz)');
    xlim([0 50]);
    
    % == first targ
    ffPBS=DATSTRUCT.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'k', 'LineStyle', '-'});
    plot(inds, ffMUSC(inds), 'sk');
    
    % == second targ
    if isfield(DATSTRUCT, 'secondtarget');
        
        ffPBS=DATSTRUCT.secondtarget(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.secondtarget(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    plot(inds, ffMUSC(inds), 'sb');
    end
    
    lt_plot_zeroline;
    
    
    
    % line for end of bidir
    if strcmp(DayToLockAllTo, 'WNday1');
        if isfield(DATSTRUCT.INFORMATION, 'BidirDay1_RelInd')
            bidir_start=DATSTRUCT.INFORMATION(i).BidirDay1_RelInd;
            if ~isempty(bidir_start)
                line([bidir_start-0.5 bidir_start-0.5], ylim, 'Color','r');
            end
        end
        bidir_last_day_actual=DATSTRUCT.INFORMATION(i).LastDay_RelInds_BidirIfExist_OtherwiseWN;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    end
    
    % line for consol start and end
    if strcmp(DayToLockAllTo, 'SingleDirConsolid');
        birdname=DATSTRUCT.INFORMATION(i).birdname;
        exptname=DATSTRUCT.INFORMATION(i).exptname;
        
        % --- need to enter consolid dates in params
        birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
        
        % convert it to relative to consol start
        consolEnd=consolPeriod(2)-consolPeriod(1)+1;
        
        line([consolEnd+0.5 consolEnd+0.5], ylim, 'Color','m');
    end
    
    
    
    %     birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    %     exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    %
    %     ind=intersect(birdInd+1, exptInd);
    %
    %     if ~isempty(ind);
    %         consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
%         
%         line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
%         line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
%     end
    
    
end

lt_subtitle('[targ + second targ] magen=consol period');


%% ====== PLOT [plus targ + mean of all same type

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT.firsttarget);
    
    birdname=DATSTRUCT.INFORMATION(i).birdname;
    exptname=DATSTRUCT.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day');
    ylabel('shift from baseline (hz)');
        xlim([0 50]);

    
    % == first targ
    ffPBS=DATSTRUCT.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'k', 'LineStyle', '-'});
    plot(inds, ffMUSC(inds), 'sk');
    
    % == second targ
    if isfield(DATSTRUCT, 'secondtarget');
        
        ffPBS=DATSTRUCT.secondtarget(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.secondtarget(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    plot(inds, ffMUSC(inds), 'sb');
    end
    
    lt_plot_zeroline;
    
    
    % ==== other nontargs (same)
    if isfield(DATSTRUCT, 'nontargets_Same');
                disp([birdname '-' exptname]);

        % === SAME
        ffPBS=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    plot(inds, ffMUSC(inds), 'sb');
    end
        
    
    
    % line for end of bidir
    if strcmp(DayToLockAllTo, 'WNday1');
        if isfield(DATSTRUCT.INFORMATION, 'BidirDay1_RelInd')
            bidir_start=DATSTRUCT.INFORMATION(i).BidirDay1_RelInd;
            if ~isempty(bidir_start)
                line([bidir_start-0.5 bidir_start-0.5], ylim, 'Color','r');
            end
        end
        bidir_last_day_actual=DATSTRUCT.INFORMATION(i).LastDay_RelInds_BidirIfExist_OtherwiseWN;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    end
    
    % line for consol start and end
    if strcmp(DayToLockAllTo, 'SingleDirConsolid');
        birdname=DATSTRUCT.INFORMATION(i).birdname;
        exptname=DATSTRUCT.INFORMATION(i).exptname;
        
        % --- need to enter consolid dates in params
        birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
        
        % convert it to relative to consol start
        consolEnd=consolPeriod(2)-consolPeriod(1)+1;
        
        line([consolEnd+0.5 consolEnd+0.5], ylim, 'Color','m');
    end
    
    
    
    %     birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    %     exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    %
    %     ind=intersect(birdInd+1, exptInd);
    %
    %     if ~isempty(ind);
    %         consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
%         
%         line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
%         line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
%     end
    
    
end

lt_subtitle('[targ + mean of same type] magen=consol period');



%% ====== PLOT [targ and other targs]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT.firsttarget);
    
    birdname=DATSTRUCT.INFORMATION(i).birdname;
    exptname=DATSTRUCT.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day');
    ylabel('shift from baseline (hz)');
    
    
    % == first targ
    ffPBS=DATSTRUCT.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'k', 'LineStyle', '-'});
%     lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Marker', 's', 'Color', 'k', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    
    % ==== other nontargs
    if isfield(DATSTRUCT, 'nontargets_Same');
        % === SAME
        ffPBS=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.nontargets_Same(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Marker', 's', 'Color', 'b', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    % === DIFF
        ffPBS=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.nontargets_Diff(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'r', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Marker', 's', 'Color', 'r', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    
    end
        
    
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    if strcmp(DayToLockAllTo, 'WNday1');
        if isfield(DATSTRUCT.INFORMATION, 'BidirDay1_RelInd')
            bidir_start=DATSTRUCT.INFORMATION(i).BidirDay1_RelInd;
            if ~isempty(bidir_start)
                line([bidir_start-0.5 bidir_start-0.5], ylim, 'Color','r');
            end
        end
        bidir_last_day_actual=DATSTRUCT.INFORMATION(i).LastDay_RelInds_BidirIfExist_OtherwiseWN;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    end
    
    % line for consol start and end
    if strcmp(DayToLockAllTo, 'SingleDirConsolid');
        birdname=DATSTRUCT.INFORMATION(i).birdname;
        exptname=DATSTRUCT.INFORMATION(i).exptname;
        
        % --- need to enter consolid dates in params
        birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
        
        % convert it to relative to consol start
        consolEnd=consolPeriod(2)-consolPeriod(1)+1;
        
        line([consolEnd+0.5 consolEnd+0.5], ylim, 'Color','m');
    end
    
    
    
    %     birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    %     exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    %
    %     ind=intersect(birdInd+1, exptInd);
    %
    %     if ~isempty(ind);
    %         consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
%         
%         line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
%         line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
%     end
    
    
end

lt_subtitle('black: targ; blue: sametype; red: difftype');



%% ====== PLOT [targ and other targs [all avg]]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(DATSTRUCT.firsttarget);
    
    birdname=DATSTRUCT.INFORMATION(i).birdname;
    exptname=DATSTRUCT.INFORMATION(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname]);
    xlabel('day');
    ylabel('shift from baseline (hz)');
    
    
    % == first targ
    ffPBS=DATSTRUCT.firsttarget(i).FFminusBase_Mean_PBS;
    ffMUSC=DATSTRUCT.firsttarget(i).FFminusBase_Mean_MUSC;
    inds=find(~isnan(ffPBS));
    
    lt_plot(inds, ffPBS(inds), {'Color', 'k', 'LineStyle', '-'});
%     lt_plot(inds, ffMUSC(inds), {'Color', 'r', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Marker', 's', 'Color', 'k', 'MarkerFaceColor','none', 'LineStyle', '-'});
    
    
    % ==== other nontargs
    if isfield(DATSTRUCT, 'nontargets_Same');
        % === ALL
        ffPBS=DATSTRUCT.nontargets_All(i).FFminusBase_Mean_PBS;
        ffMUSC=DATSTRUCT.nontargets_All(i).FFminusBase_Mean_MUSC;
        inds=find(~isnan(ffPBS));
        
    lt_plot(inds, ffPBS(inds), {'Color', 'b', 'LineStyle', '-'});
    lt_plot(inds, ffMUSC(inds), {'Marker', 's', 'Color', 'b', 'MarkerFaceColor','none', 'LineStyle', '-'});

    end
        
    
    
    
    lt_plot_zeroline;
    
    % line for end of bidir
    if strcmp(DayToLockAllTo, 'WNday1');
        if isfield(DATSTRUCT.INFORMATION, 'BidirDay1_RelInd')
            bidir_start=DATSTRUCT.INFORMATION(i).BidirDay1_RelInd;
            if ~isempty(bidir_start)
                line([bidir_start-0.5 bidir_start-0.5], ylim, 'Color','r');
            end
        end
        bidir_last_day_actual=DATSTRUCT.INFORMATION(i).LastDay_RelInds_BidirIfExist_OtherwiseWN;
        line([bidir_last_day_actual+0.5 bidir_last_day_actual+0.5], ylim, 'Color','r');
    end
    
    % line for consol start and end
    if strcmp(DayToLockAllTo, 'SingleDirConsolid');
        birdname=DATSTRUCT.INFORMATION(i).birdname;
        exptname=DATSTRUCT.INFORMATION(i).exptname;
        
        % --- need to enter consolid dates in params
        birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
        exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
        
        ind=intersect(birdInd+1, exptInd);
        
        consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};
        
        % convert it to relative to consol start
        consolEnd=consolPeriod(2)-consolPeriod(1)+1;
        
        line([consolEnd+0.5 consolEnd+0.5], ylim, 'Color','m');
    end
    
    
    
    %     birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
    %     exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
    %
    %     ind=intersect(birdInd+1, exptInd);
    %
    %     if ~isempty(ind);
    %         consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
%         
%         line([consolPeriod(1)-0.4 consolPeriod(1)-0.4], ylim, 'Color','m');
%         line([consolPeriod(2)+0.4 consolPeriod(2)+0.4], ylim, 'Color','m');
%     end
    
    
end

lt_subtitle('black: targ; blue: sametype; red: difftype');


%% ==== PLOT MEAN ACROSS ALL EXPERIMENTS FOR BIDIR
% ==== DON'T CARE ABOUT CONSOLIDATION - SIMPLY BIN IN DAYS RELATIVE TO
% LEARNING START

% === to add:
% 1) stop when reach end of bidir
% 3) only run if bird has data for all bins


maxday=length([DATSTRUCT.firsttarget.FFminusBase_Mean_PBS]);
NumDaysInBin=2;
DayBins_StartEdges=1:NumDaysInBin:maxday;

OUTPUT=struct; % one cell for each bin


for i=1:length(DayBins_StartEdges);

    OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[];
    OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[];
    OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[];
    
%     OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[];
%     OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[];
%     OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias=[];

    OUTPUT.nontargets_Same.FFrelBaseInBin(i).MUSC=[];
    OUTPUT.nontargets_Same.FFrelBaseInBin(i).PBS=[];
    OUTPUT.nontargets_Same.FFrelBaseInBin(i).MPbias=[];


    OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MUSC=[];
    OUTPUT.nontargets_Diff.FFrelBaseInBin(i).PBS=[];
    OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MPbias=[];

    OUTPUT.nontargets_All.FFrelBaseInBin(i).MUSC=[];
    OUTPUT.nontargets_All.FFrelBaseInBin(i).PBS=[];
    OUTPUT.nontargets_All.FFrelBaseInBin(i).MPbias=[];

    DaysInBin=DayBins_StartEdges(i):DayBins_StartEdges(i)+NumDaysInBin-1;

    for day=DaysInBin;

        % == for this day, collect all MUSC. for days with musc, also
        % collect PBS

        for j=1:length(DATSTRUCT.firsttarget); % num expts

            
            
            % ++++++++++++++++++++++++++++++++++++++
            % === only continue if vector is long enough to today
            if length(DATSTRUCT.firsttarget(j).FFminusBase_Mean_MUSC)<day;
                continue;
            end

            % === only continue if first target has MUSC data for today
            if isnan(DATSTRUCT.firsttarget(j).FFminusBase_Mean_MUSC(day));
                continue;
            end

            % ==== make sure this day is in consolidation period (if
            % desired)
            if strcmp(DayToLockAllTo, 'SingleDirConsolid');
            birdname=DATSTRUCT.INFORMATION(j).birdname;
            exptname=DATSTRUCT.INFORMATION(j).exptname;
                
                % --- need to enter consolid dates in params
            birdInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, birdname));
            exptInd=find(strcmp(Params.LMANTimeCourse_v2.SingleDirConsolid, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
                consolPeriod=Params.LMANTimeCourse_v2.SingleDirConsolid{ind+1};

                % convert it to relative to consol start
                consolEnd=consolPeriod(2)-consolPeriod(1)+1;
                    if day>consolEnd
                        continue
                    end
            end
            
%             if OnlyConsolPeriod_bidir==1;
%                 birdname=DATSTRUCT.INFORMATION(j).birdname;
%                 exptname=DATSTRUCT.INFORMATION(j).exptname;
% 
% 
%                 birdInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, birdname));
%                 exptInd=find(strcmp(Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart, exptname));
% 
%                 ind=intersect(birdInd+1, exptInd);
% 
%                 if ~isempty(ind);
%                     consolPeriod=Params.LMANTimeCourse.BidirConsolPeriod_IndsFromBidirStart{ind+1};
% 
%                     if day<consolPeriod(1) | day>consolPeriod(2)
%                         continue
%                     end
% 
%                 else
%                     disp(['WNARNING - no consol period entered for ' birdname '-' exptname ' - using all days']);
% 
%                 end
%             end



            % ++++++++++++++++++++++++++



            % ============== first target
            ffMUSC=DATSTRUCT.firsttarget(j).FFminusBase_Mean_MUSC(day); % MUSC
            ffPBS=DATSTRUCT.firsttarget(j).FFminusBase_Mean_PBS(day); % PBS

            % --- flip sign if learning was neg
            ffMUSC=ffMUSC*DATSTRUCT.INFORMATION(j).targ_learn_dir;
            ffPBS=ffPBS*DATSTRUCT.INFORMATION(j).targ_learn_dir;

            % ---- REVERSION (MP BIAS)
            MPbias=ffMUSC/ffPBS;

            % ---- OUTPUT
            OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC=[OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC ffMUSC];
            OUTPUT.firsttarget.FFrelBaseInBin(i).PBS=[OUTPUT.firsttarget.FFrelBaseInBin(i).PBS ffPBS];

            OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias=[OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias MPbias];

            
            % =========== OTHER TARGS -------------------------------------
            % === same type
            if ~isempty(DATSTRUCT.nontargets_Same(j).FFminusBase_Mean_MUSC)
            ffMUSC=DATSTRUCT.nontargets_Same(j).FFminusBase_Mean_MUSC(day); % MUSC
            ffPBS=DATSTRUCT.nontargets_Same(j).FFminusBase_Mean_PBS(day); % PBS

            % --- flip sign if learning was neg
            ffMUSC=ffMUSC*DATSTRUCT.INFORMATION(j).targ_learn_dir;
            ffPBS=ffPBS*DATSTRUCT.INFORMATION(j).targ_learn_dir;

            % ---- REVERSION (MP BIAS)
            MPbias=ffMUSC/ffPBS;

            % ---- OUTPUT
            OUTPUT.nontargets_Same.FFrelBaseInBin(i).MUSC=[OUTPUT.nontargets_Same.FFrelBaseInBin(i).MUSC ffMUSC];
            OUTPUT.nontargets_Same.FFrelBaseInBin(i).PBS=[OUTPUT.nontargets_Same.FFrelBaseInBin(i).PBS ffPBS];

            OUTPUT.nontargets_Same.FFrelBaseInBin(i).MPbias=[OUTPUT.nontargets_Same.FFrelBaseInBin(i).MPbias MPbias];
            end
            
            % === diff type
            if ~isempty(DATSTRUCT.nontargets_Diff(j).FFminusBase_Mean_MUSC)
            ffMUSC=DATSTRUCT.nontargets_Diff(j).FFminusBase_Mean_MUSC(day); % MUSC
            ffPBS=DATSTRUCT.nontargets_Diff(j).FFminusBase_Mean_PBS(day); % PBS

            % --- flip sign if learning was neg
            ffMUSC=ffMUSC*DATSTRUCT.INFORMATION(j).targ_learn_dir;
            ffPBS=ffPBS*DATSTRUCT.INFORMATION(j).targ_learn_dir;

            % ---- REVERSION (MP BIAS)
            MPbias=ffMUSC/ffPBS;

            % ---- OUTPUT
            OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MUSC=[OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MUSC ffMUSC];
            OUTPUT.nontargets_Diff.FFrelBaseInBin(i).PBS=[OUTPUT.nontargets_Diff.FFrelBaseInBin(i).PBS ffPBS];

            OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MPbias=[OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MPbias MPbias];
            end
            
            % === all type
            if ~isempty(DATSTRUCT.nontargets_All(j).FFminusBase_Mean_MUSC)
            ffMUSC=DATSTRUCT.nontargets_All(j).FFminusBase_Mean_MUSC(day); % MUSC
            ffPBS=DATSTRUCT.nontargets_All(j).FFminusBase_Mean_PBS(day); % PBS

            % --- flip sign if learning was neg
            ffMUSC=ffMUSC*DATSTRUCT.INFORMATION(j).targ_learn_dir;
            ffPBS=ffPBS*DATSTRUCT.INFORMATION(j).targ_learn_dir;

            % ---- REVERSION (MP BIAS)
            MPbias=ffMUSC/ffPBS;

            % ---- OUTPUT
            OUTPUT.nontargets_All.FFrelBaseInBin(i).MUSC=[OUTPUT.nontargets_All.FFrelBaseInBin(i).MUSC ffMUSC];
            OUTPUT.nontargets_All.FFrelBaseInBin(i).PBS=[OUTPUT.nontargets_All.FFrelBaseInBin(i).PBS ffPBS];

            OUTPUT.nontargets_All.FFrelBaseInBin(i).MPbias=[OUTPUT.nontargets_All.FFrelBaseInBin(i).MPbias MPbias];
            end
             % -----------------------------------------------------------------------

%             % ============= second target
%             ffMUSC=DATSTRUCT.secondtarget(j).FFminusBase_Mean_MUSC(day); % MUSC
%             ffPBS=DATSTRUCT.secondtarget(j).FFminusBase_Mean_PBS(day); % PBS
% 
%             % --- flip sign if learning was neg
%             ffMUSC=ffMUSC*DATSTRUCT.INFORMATION(j).targ_learn_dir;
%             ffPBS=ffPBS*DATSTRUCT.INFORMATION(j).targ_learn_dir;
% 
%             % ---- REVERSION (MP BIAS)
%             MPbias=ffMUSC/ffPBS;
% 
%             % ---- OUTPUT
%             OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC=[OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC ffMUSC];
%             OUTPUT.secondtarg.FFrelBaseInBin(i).PBS=[OUTPUT.secondtarg.FFrelBaseInBin(i).PBS ffPBS];
% 
%             OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias=[OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias MPbias];
        end

    end
end







% GET DATA IN DAYS RELATIVE TO CONSOLDIATION PERIOD


%% ===== PLOT (DAY BINS REL BIDIR START DAY 1)
%
% lt_figure; hold on;
% title('mean of MP biases (within each extp) (blue=second targ), (bk=first targ)');
% xlabel([ 'day bin (rel bidir day 1) (binsize: ' num2str(NumDaysInBin)]);
% ylabel('MP bias (mean of value for each expt');
%
% X=1:length(DayBins_StartEdges);
% Yfirsttarg=[];
% Ysecondtarg=[];
% Yfirsttarg_sem=[];
% Ysecondtarg_sem=[];
%
%
% for i=1:length(X);
%
%     if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias)
%     continue
%     end
%       %=== first targ
%       % -- plot each point
%       y=OUTPUT.firsttarget.FFrelBaseInBin(i).MPbias;
%
%       plot(i, y, 'ok');
%
%       % -- collect mean
%         Yfirsttarg(i)=mean(y);
%         Yfirsttarg_sem(i)=lt_sem(y);
%
%         lt_plot_bar(i, mean(y), {'Errors', lt_sem(y), 'Color','k'});
%
%         % === second targ
%        % -- plot each point
%       y=OUTPUT.secondtarg.FFrelBaseInBin(i).MPbias;
%
%       plot(i, y, 'ob');
%
%       % -- collect mean
%         Ysecondtarg(i)=mean(y);
%         Ysecondtarg_sem(i)=lt_sem(y);
%
%                 lt_plot_bar(i, mean(y), {'Errors', lt_sem(y), 'Color','b'});
%
% end
%
%
%
%% PLOT (RAW PBS AND MUSC, for targ and new targ);
lt_figure; hold on;

lt_subplot(2,1,1); hold on;
title('Raw PBS and MUSC,  all experiments');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('FF rel baseline (in targ dir)');

X=1:length(DayBins_StartEdges);
subdivisions=0.1; % for plotting diff things in one ind

for i=1:length(X);

    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end

      %=== first targ
      x=[i-subdivisions*5 i-subdivisions*4];
      % -- plot each point
      yPBS=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
      yMUSC=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;

      plot(x, [yPBS; yMUSC], '-ok');

      % -- means
      if size([yPBS' yMUSC'],1)>1
%       lt_plot_bar(x-0.05, mean([yPBS' yMUSC']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])});
      lt_plot_bar(x-subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'k'});
      end

      % ====== NON TARGS --------------------
      % == same type
      x=[i-subdivisions*2 i-subdivisions];
      % -- plot each point
      yPBS=OUTPUT.nontargets_Same.FFrelBaseInBin(i).PBS;
      yMUSC=OUTPUT.nontargets_Same.FFrelBaseInBin(i).MUSC;

%       plot(x, [yPBS' yMUSC'], '-ob');
      plot(x, [yPBS; yMUSC], '-ob');

      % -- means
      if size([yPBS' yMUSC'],1)>1
%       lt_plot_bar(x-0.05, mean([yPBS' yMUSC']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])});
      lt_plot_bar(x-subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'b'});
      end
      
      
      % == diff type
      x=[i i+subdivisions];
      % -- plot each point
      yPBS=OUTPUT.nontargets_Diff.FFrelBaseInBin(i).PBS;
      yMUSC=OUTPUT.nontargets_Diff.FFrelBaseInBin(i).MUSC;

%       plot(x, [yPBS' yMUSC'], '-ob');
      plot(x, [yPBS; yMUSC], '-or');

      % -- means
      if size([yPBS' yMUSC'],1)>1
%       lt_plot_bar(x-0.05, mean([yPBS' yMUSC']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])});
      lt_plot_bar(x-subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'r'});
      end
      
        % === all
       x=[i+subdivisions*2 i+subdivisions*3];
      % -- plot each point
      yPBS=OUTPUT.nontargets_All.FFrelBaseInBin(i).PBS;
      yMUSC=OUTPUT.nontargets_All.FFrelBaseInBin(i).MUSC;

%       plot(x, [yPBS' yMUSC'], '-ob');
      plot(x, [yPBS; yMUSC], '-om');

      % -- means
      if size([yPBS' yMUSC'],1)>1
%       lt_plot_bar(x-0.05, mean([yPBS' yMUSC']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])});
      lt_plot_bar(x-subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'm'});
      end
       
        lt_plot_annotation('', {'bk=targ','bl=same','red=diff','mag=all'},'k');
      
      
%       % === second targ
%       x=[i+subdivisions i+2*subdivisions];
% 
%       % -- plot each point
%       yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
%       yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
% 
%       plot(x, [yPBS' yMUSC'], '-ob');
% 
%       % -- means
%       if size([yPBS' yMUSC'],1)>1
% %       lt_plot_bar(x+0.05, mean([yPBS' yMUSC']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])});
%       lt_plot_bar(x+subdivisions/4, mean([yPBS' yMUSC']), {'Color', 'b'});
%       end

end


% ========= [ only mean, each mean norm to PBS]
lt_subplot(2,1,2); hold on;
title('Raw PBS and MUSC, all experiments [mean, norm so PBS is 1]');
xlabel([ 'day bin (binsize: ' num2str(NumDaysInBin) ')']);
ylabel('FF rel baseline (in targ dir)');

for i=1:length(X);

    if isempty(OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC)
        continue
    end

      %=== first targ
      x=[i-subdivisions*2 i-subdivisions];
      % -- plot each point
      yPBS=OUTPUT.firsttarget.FFrelBaseInBin(i).PBS;
      yMUSC=OUTPUT.firsttarget.FFrelBaseInBin(i).MUSC;


      % -- means
      if size([yPBS' yMUSC'],1)>1
      lt_plot_bar(x-0.05, mean([yPBS' yMUSC'])./mean([yPBS' yPBS']), {'Color', 'k', 'Errors', lt_sem([yPBS' yMUSC'])./mean([yPBS' yPBS'])});
      end


%       % === second targ
%       x=[i+subdivisions i+2*subdivisions];
% 
%       % -- plot each point
%       yPBS=OUTPUT.secondtarg.FFrelBaseInBin(i).PBS;
%       yMUSC=OUTPUT.secondtarg.FFrelBaseInBin(i).MUSC;
% 
%       % -- means
%       if size([yPBS' yMUSC'],1)>1
%       lt_plot_bar(x-0.05, -mean([yPBS' yMUSC'])./mean([yPBS' yPBS']), {'Color', 'b', 'Errors', lt_sem([yPBS' yMUSC'])./mean([yPBS' yPBS'])});
%       end

end





% %% ==================== COLLECT, BUT HAVE FIRST WN DAY BE LOCKED TO DAY 1 (INCLUDE BIDIR, BUT ANNOTATE AS SUCH)
%


