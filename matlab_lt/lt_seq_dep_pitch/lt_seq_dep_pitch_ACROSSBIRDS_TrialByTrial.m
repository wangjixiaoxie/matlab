%% LT 1/22/16 -
function lt_seq_dep_pitch_ACROSSBIRDS_TrialByTrial(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, OnlyWellLearned)

if OnlyWellLearned==1;
    filter='learning_metric';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
end
    
%% ==

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% ++++++++++++++++++++ PLOTS OVER ALL BIRDS
        figcount=1;
        subplotrows=4;
        subplotcols=3;
        fignums_alreadyused=[];
        hfigs=[];
        hsplot_all=[];
        
%% ===== PLOT RELATIONSHIP BETWEEN CORR OF WITHIN DAY SLOPES OVER LEARNING AND GENERALIZATION

Generalization_AllBirds=[];
DayToDayCorrOfSlopes_AllBirds=[];
Similar_AllBirds=[];
MeanSlopeAdaptationTargAllBirds=[];
MeanSlopeAdaptationAllBirds=[];
         MeanSlopeConsolTargAllBirds=[];
        MeanSlopeConsolAllBirds=[];


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            datafield_ff='FFvals_DevFromBase_WithinTimeWindow';
            datafield_tvals='Tvals_WithinTimeWindow';
        else
            datafield_ff='FFvals';
            datafield_tvals='Tvals';
        end
        
        
        
        % ====== COLLECT DAY TO DAY SLOPE FOR TARGET ======================
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        syl=targsyl;
        numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
            {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
        
        SlopeAllDaysThisSyl=[];
        for k=1:numdays
            
            if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                continue
            end
            
            % === for this day, collect other stats (e.g. slope)
            ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
            tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
            
            if iscell(ffvals)
                ffvals=cell2mat(ffvals);
            end
            if iscell(tvals)
                tvals=cell2mat(tvals);
            end
            
            % convert tvals to hours within this day
            [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
            tvals=tmp2.hours;
            
            [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
            
            % ===== OUTPUT
            
            SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
            
            
            %                 DATSTRUCT.experiment(ind).syllable(ind)
        end
%         SlopeAllDaysThisSyl_Targ=SlopeAllDaysThisSyl_Targ.*SeqDepPitch_AcrossBirds.birds{i}.experiment...
%             {ii}.INFORMATION.targ_learn_dir;
        SlopeAllDaysThisSyl_Targ=SlopeAllDaysThisSyl;
        
        
        
        % ================== COLLECT FOR ALL OTHER SYLS
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        GeneralizationAll=[];
        DayToDayCorrOfSlopes=[];
        DayToDayCorrOfSlopes_base=[];
        SimilarAll=[];
        PreSimAll=[];
        
        MeanSlopeAdaptationAll=[];
        MeanSlopeAdaptationTargAll=[];
                        MeanSlopeConsolAll=[];
                MeanSlopeConsolTargAll=[];

        
        % ===== Collect trial by trail pitch and times for each
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            %             if istarget
            %                 continue;
            %             end
            
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
            
            
            SlopeAllDaysThisSyl=[];
            for k=1:numdays
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                    SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                    continue
                end
                
                % === for this day, collect other stats (e.g. slope)
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
                
                if iscell(ffvals)
                    ffvals=cell2mat(ffvals);
                end
                if iscell(tvals)
                    tvals=cell2mat(tvals);
                end
                
                % convert tvals to hours within this day
                [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
                tvals=tmp2.hours;
                
                [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
                
                % ===== OUTPUT
                
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
                
            end
            
%             SlopeAllDaysThisSyl=SlopeAllDaysThisSyl.*SeqDepPitch_AcrossBirds.birds{i}.experiment...
%                 {ii}.INFORMATION.targ_learn_dir;
%             
            
            % ======= COLLECT MEAN OVER BASELINE, ADAPTIVE, AND CONSOLID
            % DAYS
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            BaselineDays=1:WNday1-1;
            SlopeMeanBaseline=nanmean(SlopeAllDaysThisSyl(BaselineDays));
            
            % get consolidation days
            birdInd=find(strcmp(PARAMS.TimeCourse.ConsolidationList, birdname));
            exptInd=find(strcmp(PARAMS.TimeCourse.ConsolidationList, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
            if ~isempty(ind);
                consolPeriod=PARAMS.TimeCourse.ConsolidationList{ind+1};
            else
                disp('PROBLEM - no consolid inds for this bird');
                lt_plot_text(1, 1, 'PROBLEM - no consolid inds for this bird');
                continue
            end
            
            
            consolPeriod_OverallInds=consolPeriod+WNday1-1;
            
            AdaptiveDays=WNday1:consolPeriod_OverallInds-1;
            ConsolidDays=consolPeriod_OverallInds(1):consolPeriod_OverallInds(2);
            WNDays=AdaptiveDays(1):ConsolidDays(end);
            
            SlopeMeanAdaptive=nanmean(SlopeAllDaysThisSyl(AdaptiveDays));
            SlopeMeanConsolid=nanmean(SlopeAllDaysThisSyl(ConsolidDays));
            
            SlopeMeanWNdays=nanmean(SlopeAllDaysThisSyl(WNDays));
            
            
            % ==== what is day to day correlation in slopes to the target?
            WNDays2=WNDays(~isnan(SlopeAllDaysThisSyl(WNDays)));
            SlopeCorrVsTarg=corr(SlopeAllDaysThisSyl(WNDays2)', SlopeAllDaysThisSyl_Targ(WNDays2)');
            
            BaselineDays2=BaselineDays(~isnan(SlopeAllDaysThisSyl(BaselineDays)));
            SlopeCorrVsTarg_base=corr(SlopeAllDaysThisSyl(BaselineDays2)', SlopeAllDaysThisSyl_Targ(BaselineDays2)');
            
            % === what is learning genearlization?
            GeneralizationLearning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            
            
            
            
            % ===== OUTPUT for 45 degree scatter
            if ~istarget
                SimilarAll=[SimilarAll similar];
                PreSimAll=[PreSimAll presimilar];
                GeneralizationAll=[GeneralizationAll GeneralizationLearning];
                DayToDayCorrOfSlopes=[DayToDayCorrOfSlopes SlopeCorrVsTarg];
                DayToDayCorrOfSlopes_base=[DayToDayCorrOfSlopes_base SlopeCorrVsTarg_base];
                
                MeanSlopeAdaptationAll=[MeanSlopeAdaptationAll SlopeMeanAdaptive];
                MeanSlopeAdaptationTargAll=[MeanSlopeAdaptationTargAll nanmean(SlopeAllDaysThisSyl_Targ(AdaptiveDays))];

                MeanSlopeConsolAll=[MeanSlopeConsolAll SlopeMeanConsolid];
                MeanSlopeConsolTargAll=[MeanSlopeConsolTargAll nanmean(SlopeAllDaysThisSyl_Targ(ConsolidDays))];

            end
            
        end
        
        
        % ===== PLOT FOR THISE XPT
        if (0)
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplot_all=[hsplot_all hsplot];
            title([birdname '-' exptname]);
            
            ylabel('corr to target, day to day slopes (WN)');
            xlabel('generalization');
            
            % same type
            color='b';
            inds=SimilarAll==1;
            X=GeneralizationAll(inds);
            Y=DayToDayCorrOfSlopes(inds);
            if ~isempty(X)
                lt_plot_45degScatter(X, Y, color);
            end
            
            % diff type
            color='r';
            inds=SimilarAll==0;
            X=GeneralizationAll(inds);
            Y=DayToDayCorrOfSlopes(inds);
            
            if ~isempty(X)
                lt_plot_45degScatter(X, Y, color);
            end
        end
        
        
        % ===== COLLECT ACROSS ALL BIRDS
        Generalization_AllBirds=[Generalization_AllBirds GeneralizationAll];
        DayToDayCorrOfSlopes_AllBirds=[DayToDayCorrOfSlopes_AllBirds DayToDayCorrOfSlopes];
        Similar_AllBirds=[Similar_AllBirds SimilarAll];
        
        MeanSlopeAdaptationTargAllBirds=[MeanSlopeAdaptationTargAllBirds MeanSlopeAdaptationTargAll];
        MeanSlopeAdaptationAllBirds=[MeanSlopeAdaptationAllBirds MeanSlopeAdaptationAll];
 
         MeanSlopeConsolTargAllBirds=[MeanSlopeConsolTargAllBirds MeanSlopeConsolTargAll];
        MeanSlopeConsolAllBirds=[MeanSlopeConsolAllBirds MeanSlopeConsolAll];

    end
end



%% ==================== PLOT ACROSS ALL EXPERIMENTS - PLOT 1
lt_figure; hold on;
        % ===== PLOT FOR THISE XPT
        ylabel('corr to target, day to day slopes (WN)');
        xlabel('generalization');
        
        % same type
        color='b';
        inds=Similar_AllBirds==1;
        X=Generalization_AllBirds(inds);
        Y=DayToDayCorrOfSlopes_AllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
        % diff type
        color='r';
        inds=Similar_AllBirds==0;
        X=Generalization_AllBirds(inds);
        Y=DayToDayCorrOfSlopes_AllBirds(inds);
        
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end


% ==================== PLOT ACROSS ALL EXPERIMENTS - PLOT 2 - SLOPE DURING
% ADAPTATION
lt_figure; hold on;
        % ===== PLOT FOR THISE XPT
        xlabel('target mean slope (adaptation period)');
        ylabel('nontarget mean slope (adaptation period)');
        
        % same type
        color='b';
        inds=Similar_AllBirds==1;
        Y=MeanSlopeAdaptationAllBirds(inds);
        X=MeanSlopeAdaptationTargAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
        % diff type
        color='r';
        inds=Similar_AllBirds==0;
        Y=MeanSlopeAdaptationAllBirds(inds);
        X=MeanSlopeAdaptationTargAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
    % ==================== PLOT ACROSS ALL EXPERIMENTS - PLOT 2 - SLOPE DURING
% CONSOL
lt_figure; hold on;
        % ===== PLOT FOR THISE XPT
        xlabel('target mean slope (consol period)');
        ylabel('nontarget mean slope (consol period)');
        
        % same type
        color='b';
        inds=Similar_AllBirds==1;
        Y=MeanSlopeConsolAllBirds(inds);
        X=MeanSlopeConsolTargAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
        % diff type
        color='r';
        inds=Similar_AllBirds==0;
        Y=MeanSlopeConsolAllBirds(inds);
        X=MeanSlopeConsolTargAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
      % ==================== PLOT ACROSS ALL EXPERIMENTS - PLOT 2 - SLOPE DURING
% ADAPT VS CONSOL
lt_figure; hold on;
        % ===== PLOT FOR THISE XPT
        xlabel('nontarg slope during consol');
        ylabel('nontarg slope during adapt');
        
        % same type
        color='b';
        inds=Similar_AllBirds==1;
        Y=MeanSlopeAdaptationAllBirds(inds);
        X=MeanSlopeConsolAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
        
        % diff type
        color='r';
        inds=Similar_AllBirds==0;
        Y=MeanSlopeAdaptationAllBirds(inds);
        X=MeanSlopeConsolAllBirds(inds);
        if ~isempty(X)
            lt_plot_45degScatter(X, Y, color);
        end
  
        
pause;




%% ===== plot 45 degree  scatter (generalization vs day to day corr of slopes)
lt_figure; hold on;     
ylabel('corr to target, day to day slopes (WN)');
xlabel('generalization');

% same type
color='b';
inds=SimilarAll==1;
X=GeneralizationAll(inds);
Y=DayToDayCorrOfSlopes(inds);
if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end

% diff type
color='r';
inds=SimilarAll==0;
X=GeneralizationAll(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end


% ====================
lt_figure; hold on;     
ylabel('corr to target, day to day slopes (WN)');
xlabel('corr to target, day to day slopes (base)');

% same type
color='b';
inds=SimilarAll==1;
X=DayToDayCorrOfSlopes_base(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end
            
% diff type
color='r';
inds=SimilarAll==0;
X=DayToDayCorrOfSlopes_base(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end



%% ++++++++++++++++++++++++++++ ASINGLE BIRD STUFF
%% ===== PLOT MEANS OF SLOPES DURING BASELINE, ADAPTIVE, AND CONSOLID PHASES
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if ~strcmp(birdname, BirdToPlot) | ~strcmp(exptname, ExptToPlot);
            continue
        end
        
        figcount=1;
        subplotrows=4;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplot_all=[];
        
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            datafield_ff='FFvals_DevFromBase_WithinTimeWindow';
            datafield_tvals='Tvals_WithinTimeWindow';
        else
            datafield_ff='FFvals';
            datafield_tvals='Tvals';
        end
        
        % ====== COLLECT DAY TO DAY SLOPE FOR TARGET ======================
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            syl=targsyl;
            numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
                        
            SlopeAllDaysThisSyl=[];
            for k=1:numdays
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                    SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                    continue
                end
                
                % === for this day, collect other stats (e.g. slope)
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
                
                if iscell(ffvals)
                    ffvals=cell2mat(ffvals);
                end
                if iscell(tvals)
                    tvals=cell2mat(tvals);
                end
                
                % convert tvals to hours within this day
                [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
                tvals=tmp2.hours;
                
                [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
                
                % ===== OUTPUT
                
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
                
                
                %                 DATSTRUCT.experiment(ind).syllable(ind)
            end
            SlopeAllDaysThisSyl_Targ=SlopeAllDaysThisSyl;
            
            
            
        % ================== COLLECT FOR ALL OTHER SYLS
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        GeneralizationAll=[];
        DayToDayCorrOfSlopes=[];
        DayToDayCorrOfSlopes_base=[];
        SimilarAll=[];
        PreSimAll=[];
        
        % ===== Collect trial by trail pitch and times for each
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
%             if istarget
%                 continue;
%             end
%             
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
            
            
            SlopeAllDaysThisSyl=[];
            for k=1:numdays
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                    SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                    continue
                end
                
                % === for this day, collect other stats (e.g. slope)
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
                
                if iscell(ffvals)
                    ffvals=cell2mat(ffvals);
                end
                if iscell(tvals)
                    tvals=cell2mat(tvals);
                end
                
                % convert tvals to hours within this day
                [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
                tvals=tmp2.hours;
                
                [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
                
                % ===== OUTPUT
                
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
                
            end
            
            
            
            % ======= COLLECT MEAN OVER BASELINE, ADAPTIVE, AND CONSOLID
            % DAYS
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            BaselineDays=1:WNday1-1;
            SlopeMeanBaseline=nanmean(SlopeAllDaysThisSyl(BaselineDays));
            
            % get consolidation days
            birdInd=find(strcmp(PARAMS.TimeCourse.ConsolidationList, birdname));
            exptInd=find(strcmp(PARAMS.TimeCourse.ConsolidationList, exptname));
            
            ind=intersect(birdInd+1, exptInd);
            
            if ~isempty(ind);
                consolPeriod=PARAMS.TimeCourse.ConsolidationList{ind+1};
            else
                disp('PROBLEM - no consolid inds for this bird');
                eadcjaeiownec;
            end
            
            
            consolPeriod_OverallInds=consolPeriod+WNday1-1;
           
            AdaptiveDays=WNday1:consolPeriod_OverallInds-1;
            ConsolidDays=consolPeriod_OverallInds(1):consolPeriod_OverallInds(2);
            WNDays=AdaptiveDays(1):ConsolidDays(2);
            
            SlopeMeanAdaptive=nanmean(SlopeAllDaysThisSyl(AdaptiveDays));
            SlopeMeanConsolid=nanmean(SlopeAllDaysThisSyl(ConsolidDays));
            
            SlopeMeanWNdays=nanmean(SlopeAllDaysThisSyl(WNDays));
            
             
            % ==== what is day to day correlation in slopes to the target?
            WNDays2=WNDays(~isnan(SlopeAllDaysThisSyl(WNDays)));
            SlopeCorrVsTarg=corr(SlopeAllDaysThisSyl(WNDays2)', SlopeAllDaysThisSyl_Targ(WNDays2)');
            
            BaselineDays2=BaselineDays(~isnan(SlopeAllDaysThisSyl(BaselineDays)));
            SlopeCorrVsTarg_base=corr(SlopeAllDaysThisSyl(BaselineDays2)', SlopeAllDaysThisSyl_Targ(BaselineDays2)');
            
            % === what is learning genearlization?
            GeneralizationLearning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            
            
            
            % ==== PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplot_all=[hsplot_all hsplot];
            title(syl);
            
            if istarget==1;
                color='k';
            elseif similar & presimilar
                color='b';
            elseif similar & ~presimilar
                color='c';
            elseif ~similar & presimilar
                color='r';
            elseif ~similar & ~presimilar
                color='m';
            end
            
            % ===== bar plot
            X=1:5;
            Y=[SlopeMeanBaseline SlopeMeanAdaptive SlopeMeanConsolid SlopeMeanWNdays SlopeCorrVsTarg];
            
            lt_plot_bar(X, Y, {'Color',color});
            
            % ===== OUTPUT for 45 degree scatter
            if ~istarget
                SimilarAll=[SimilarAll similar];
                PreSimAll=[PreSimAll presimilar];
            GeneralizationAll=[GeneralizationAll GeneralizationLearning];
            DayToDayCorrOfSlopes=[DayToDayCorrOfSlopes SlopeCorrVsTarg];
            DayToDayCorrOfSlopes_base=[DayToDayCorrOfSlopes_base SlopeCorrVsTarg_base];
            
            end
                
        end        
        
        % =======
        linkaxes(hsplot_all,'xy');
        
    end
end

lt_subtitle('slope(base), slope(adapt), slope(consol), slope (WN), corrTarg');


%% ===== plot 45 degree  scatter (generalization vs day to day corr of slopes)
lt_figure; hold on;     
ylabel('corr to target, day to day slopes (WN)');
xlabel('generalization');

% same type
color='b';
inds=SimilarAll==1;
X=GeneralizationAll(inds);
Y=DayToDayCorrOfSlopes(inds);
if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end

% diff type
color='r';
inds=SimilarAll==0;
X=GeneralizationAll(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end


% ====================
lt_figure; hold on;     
ylabel('corr to target, day to day slopes (WN)');
xlabel('corr to target, day to day slopes (base)');

% same type
color='b';
inds=SimilarAll==1;
X=DayToDayCorrOfSlopes_base(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end
            
% diff type
color='r';
inds=SimilarAll==0;
X=DayToDayCorrOfSlopes_base(inds);
Y=DayToDayCorrOfSlopes(inds);

if ~isempty(X)
lt_plot_45degScatter(X, Y, color);
end



%% === FOR EACH SYL PLOT SLOPE DURING DAY
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if ~strcmp(birdname, BirdToPlot) | ~strcmp(exptname, ExptToPlot);
            continue
        end
        
        figcount=1;
        subplotrows=4;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplot_all=[];
        
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            datafield_ff='FFvals_DevFromBase_WithinTimeWindow';
            datafield_tvals='Tvals_WithinTimeWindow';
        else
            datafield_ff='FFvals';
            datafield_tvals='Tvals';
        end
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        
        % ===== Collect trial by trail pitch and times for each
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
            
            
            SlopeAllDaysThisSyl=[];
            for k=1:numdays
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                    SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                    continue
                end
                
                % === for this day, collect other stats (e.g. slope)
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
                
                if iscell(ffvals)
                    ffvals=cell2mat(ffvals);
                end
                if iscell(tvals)
                    tvals=cell2mat(tvals);
                end
                
                % convert tvals to hours within this day
                [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
                tvals=tmp2.hours;
                
                [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
                
                % ===== OUTPUT
                
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
                
                
                %                 DATSTRUCT.experiment(ind).syllable(ind)
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplot_all=[hsplot_all hsplot];
            title(syl);
            
            if istarget==1;
                color='k';
            elseif similar & presimilar
                color='b';
            elseif similar & ~presimilar
                color='c';
            elseif ~similar & presimilar
                color='r';
            elseif ~similar & ~presimilar
                color='m';
            end
            
            
            lt_plot_bar(1:length(SlopeAllDaysThisSyl), SlopeAllDaysThisSyl, {'Color',color});
            
            
            % ============== annotate with lines for stuff
            Ylim=ylim;
            Yrange=Ylim(2)-Ylim(1);
            Xlim=xlim;
            
            % WN on and off
            WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
            WNoffday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            line([WNoffday+0.5 WNoffday+0.5],ylim,'Color','r');
            lt_plot_text(WNonday, Ylim(1), 'WNon/off','r')
            
            % day I called consol start
            try
                day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd;
                
                line([day1-0.5 day1-0.5], ylim);
                line([day2+0.5 day2+0.5], ylim);
                
                lt_plot_text(day1, Ylim(1), 'consol','b')
            catch err
            end
            
            % multidir
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
                date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
                date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
                
                line([date1_ind+0.4 date1_ind+0.4], ylim, 'Color', 'k');
                line([date2_ind+0.4 date2_ind+0.4], ylim, 'Color', 'k');
                lt_plot_text(date1_ind, Ylim(1), 'multidir','k')
            end
            
            % LMAN inactivation days
            if  SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
                targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
                
                % plot
                plot(muscdays_inds, Ylim(1)+Yrange/6, '^k', 'MarkerSize', 8);
                plot(muscdays_inds, Ylim(2)-Yrange/6, 'vk', 'MarkerSize', 8);
                
                %             for k=1:length(muscdays_inds);
                % %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'Color',[0.7 0.7 0.7]);
                %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'LineStyle','--','Color',[0.3 0.3 0.3]);
                %             end
                
                lt_plot_text(1, Ylim(1)+Yrange/8, 'arrowhead: MUSC', 'k');
                
                
            end
            
            % annotate target learning direction (as detected automatically)
            targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            lt_plot_text(2, Ylim(2)-Yrange/8, ['targ learn dir: ' num2str(targ_learn_dir)], 'r');
            
            
            
        end
        
        
        
        
        
        % =======
        linkaxes(hsplot_all,'xy');
        
    end
end

%% === FOR EACH SYL PLOT O/N CHANGE (USING 20 RENDITIONS PRE AND POST)
if (0)
Nrends=20;


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if ~strcmp(birdname, BirdToPlot) | ~strcmp(exptname, ExptToPlot);
            continue
        end
        
        figcount=1;
        subplotrows=4;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplot_all=[];
        
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            datafield_ff='FFvals_DevFromBase_WithinTimeWindow';
            datafield_tvals='Tvals_WithinTimeWindow';
        else
            datafield_ff='FFvals';
            datafield_tvals='Tvals';
        end
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        
        % ===== Collect trial by trail pitch and times for each
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            numdays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff));
            
            
            SlopeAllDaysThisSyl=[];
            for k=1:numdays
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k})
                    SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl nan];
                    continue
                end
                
                % === for this day, collect other stats (e.g. slope)
                ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){k};
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){k};
                
                if iscell(ffvals)
                    ffvals=cell2mat(ffvals);
                end
                if iscell(tvals)
                    tvals=cell2mat(tvals);
                end
                
                % convert tvals to hours within this day
                [ ~, tmp2]=lt_convert_datenum_to_hour(tvals);
                tvals=tmp2.hours;
                
                [~, ~, ~, ~, ~, SummaryStats]=lt_regress(ffvals, tvals, 0);
                
                % ===== OUTPUT
                
                SlopeAllDaysThisSyl=[SlopeAllDaysThisSyl SummaryStats.slope];
                
                
                %                 DATSTRUCT.experiment(ind).syllable(ind)
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplot_all=[hsplot_all hsplot];
            title(syl);
            
            if istarget==1;
                color='k';
            elseif similar & presimilar
                color='b';
            elseif similar & ~presimilar
                color='c';
            elseif ~similar & presimilar
                color='r';
            elseif ~similar & ~presimilar
                color='m';
            end
            
            
            lt_plot_bar(1:length(SlopeAllDaysThisSyl), SlopeAllDaysThisSyl, {'Color',color});
            
            
            % ============== annotate with lines for stuff
            Ylim=ylim;
            Yrange=Ylim(2)-Ylim(1);
            Xlim=xlim;
            
            % WN on and off
            WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
            WNoffday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            line([WNoffday+0.5 WNoffday+0.5],ylim,'Color','r');
            lt_plot_text(WNonday, Ylim(1), 'WNon/off','r')
            
            % day I called consol start
            try
                day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd;
                
                line([day1-0.5 day1-0.5], ylim);
                line([day2+0.5 day2+0.5], ylim);
                
                lt_plot_text(day1, Ylim(1), 'consol','b')
            catch err
            end
            
            % multidir
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
                date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
                date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
                
                line([date1_ind+0.4 date1_ind+0.4], ylim, 'Color', 'k');
                line([date2_ind+0.4 date2_ind+0.4], ylim, 'Color', 'k');
                lt_plot_text(date1_ind, Ylim(1), 'multidir','k')
            end
            
            % LMAN inactivation days
            if  SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
                targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
                
                % plot
                plot(muscdays_inds, Ylim(1)+Yrange/6, '^k', 'MarkerSize', 8);
                plot(muscdays_inds, Ylim(2)-Yrange/6, 'vk', 'MarkerSize', 8);
                
                %             for k=1:length(muscdays_inds);
                % %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'Color',[0.7 0.7 0.7]);
                %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'LineStyle','--','Color',[0.3 0.3 0.3]);
                %             end
                
                lt_plot_text(1, Ylim(1)+Yrange/8, 'arrowhead: MUSC', 'k');
                
                
            end
            
            % annotate target learning direction (as detected automatically)
            targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            lt_plot_text(2, Ylim(2)-Yrange/8, ['targ learn dir: ' num2str(targ_learn_dir)], 'r');
            
            
        end
        
        
        
        
        
        % =======
        linkaxes(hsplot_all,'xy');
        
    end
end




end
end
