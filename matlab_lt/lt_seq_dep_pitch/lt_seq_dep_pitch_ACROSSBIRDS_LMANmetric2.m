function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANbidir(SeqDepPitch_AcrossBirds, PARAMS, use_final_extracted_windows, NormToTarg)
%% LT 1/31/16 - BIDIR LMAN, CALCULATING A SEPARATION METRIC (AFP INCREASE SEPARATION MORE THAN MP?)


%% SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA

filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

TotalNumExperiments=0;
for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
TotalNumExperiments=TotalNumExperiments+numexperiments;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_BiDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_BiDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_BiDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    if ~any(ind_to_remove)
                        % then syl not found...
                        disp('Problem - syl not found for inds to remove');
                        continue; 
                    end
                    
                    syls_unique(ind_to_remove)=[];
                    
                    % Put back into structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=syls_unique;
                    
                    % tell user this is done
                    disp([birdname '-' exptname ': removed ' tmp_sylremove]);
                end
            end
        end
    end
end



%% ===== [COLLECT] PLOT ALL IN ONE PLOT
figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

Index1_All=[];
SameType_All=[];

DATSTRUCT.targ1.SeparationLearningEarly=[];
DATSTRUCT.targ1.SeparationMPearly=[];
DATSTRUCT.targ1.SeparationAFPearly=[];

DATSTRUCT.targ1.SeparationLearningLate=[];
DATSTRUCT.targ1.SeparationMPlate=[];
DATSTRUCT.targ1.SeparationAFPlate=[];

DATSTRUCT.targ2.SeparationLearningEarly_All=[];
DATSTRUCT.targ2.SeparationMPearly=[];
DATSTRUCT.targ2.SeparationAFPearly=[];

DATSTRUCT.targ2.SeparationLearningLate=[];
DATSTRUCT.targ2.SeparationMPlate=[];
DATSTRUCT.targ2.SeparationAFPlate=[];


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== MAKE SURE THIS HAS BIDIR LEARNING
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'days_bidir_early') & ~isfield(...
                SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'final_extracted_window_bidir')
            continue
        end
        
        
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        targsyl_other=tmp{~strcmp(tmp,targsyl)}; % find the other syl assuming that multidirsyls contains the target and the other target
        rng('shuffle');
        
        % == Only continue if both targsyl and targsyl_other are part of
        % SylsUnique
        if ~any(strcmp(SylsUnique, targsyl)) | ~any(strcmp(SylsUnique, targsyl_other));
            disp(['Skipping ' birdname '-' exptname ' - one or other targsyl not part of SylsUnique']);
            continue
        end
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_consolid_late');
            window1='days_consolid_late';
            
        elseif isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_consolid_early')
            window1 = 'days_consolid_early';            
        else
            window1='final_extracted_window';
            disp('USING EXTRACTED WINDOW FOR EARLY !!! NOT GOOD!!');
        end
        
        if use_final_extracted_windows==1;
            EpochNameList={window1, 'final_extracted_window_bidir'};
        else
            % USE BIDIR LATE. If don't have, then use bidir early
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'days_bidir_late');
                EpochNameList={window1, 'days_bidir_late'};
            else
                EpochNameList={window1, 'days_bidir_early'};
            end
        end
        
        % ============ CALCULATE METRICS
        % 1) Metric: Index 1, asks how much of the increase in separation
        % at target learning is accounted for by increase in AFP separation
        % vs. increase in MP separation. null is 0, positive is more AFP
        
        % --- early
        epochfield=EpochNameList{1};
        
        FF_PBS_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs; % mean using rends across days
        FF_MUSC_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
        FF_AFP_targ1=FF_PBS_targ1-FF_MUSC_targ1;
        
        
        FF_PBS_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl_other).meanFF_pbs; % mean using rends across days
        FF_MUSC_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl_other).meanFF_musc;
        FF_AFP_targ2=FF_PBS_targ2-FF_MUSC_targ2;
        
        SeparationLearningEarly=FF_PBS_targ1-FF_PBS_targ2;
        SeparationMPearly=FF_MUSC_targ1-FF_MUSC_targ2;
        SeparationAFPearly=FF_AFP_targ1-FF_AFP_targ2;
        
         % --- late
        epochfield=EpochNameList{2};
        
        FF_PBS_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs; % mean using rends across days
        FF_MUSC_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
        FF_AFP_targ1=FF_PBS_targ1-FF_MUSC_targ1;
        
        
        FF_PBS_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl_other).meanFF_pbs; % mean using rends across days
        FF_MUSC_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl_other).meanFF_musc;
        FF_AFP_targ2=FF_PBS_targ2-FF_MUSC_targ2;
        
        SeparationLearningLate=FF_PBS_targ1-FF_PBS_targ2;
        SeparationMPlate=FF_MUSC_targ1-FF_MUSC_targ2;
        SeparationAFPlate=FF_AFP_targ1-FF_AFP_targ2;
        
        
        % ---- CALCULATE
        % if this is >0, then means MP increased separation
        indexMP=(SeparationMPlate/SeparationLearningLate) ...
            - (SeparationMPlate + SeparationMPearly)/(SeparationLearningLate+SeparationLearningEarly);
        
        indexAFP=(SeparationAFPlate/SeparationLearningLate) ...
            - (SeparationAFPlate + SeparationAFPearly)/(SeparationLearningLate+SeparationLearningEarly);
        
        Index1=indexAFP-indexMP; % if this >0, then means AFP increased separation more than MP increased separation
       
        
        % ======================= COLLECT
        Index1_All=[Index1_All Index1];
        
        sametype=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl_other).similar_to_targ;
        SameType_All=[SameType_All sametype];
        
        
        
        % ===== PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['index1: ' num2str(Index1)]);
        lt_plot_annotation(1, [birdname '-' exptname '-(targ) ' targsyl '-(other) ' targsyl_other],'k');
        
        Y=[SeparationLearningEarly, SeparationLearningLate, SeparationMPearly, ...
            SeparationMPlate, SeparationAFPearly, SeparationAFPlate];
        X=[1 2 4 5 7 8];
        lt_plot_bar(X, Y)

        
        
    end
end

figure(hfigs);
lt_plot_annotation(1, 'x bars: early vs. late; in order: Learn, MP, AFP', 'k');


            
            
            %% +========= PLOT
            
            lt_figure; hold on;
            title('same type, index 1');
            
            inds=SameType_All==1;
            plot(1, Index1_All(inds), 'ok');
            line(xlim, [0 0]);

            
             lt_figure; hold on;
            title('all type, index 1');

            plot(1, Index1_All, 'ok');
            line(xlim, [0 0]);

           
            
            
            
            