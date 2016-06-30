function    lt_seq_dep_pitch_ACROSSBIRDS_LMANseparation(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl);
%% LT 10/12/15 - run this from _LMANlearning


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

disp('--');
disp('Removing syllables that should not be analyzed - i.e. WN overlap, since not using catch songs. REMOVED:');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_SingleDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
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


%% FOR EACH EXPERIMENT, plot AFP inactivation for each day

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
    % INITIATE FIGS
    count=1;
    SubplotsPerFig=16;
    subplotrows=4;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        DaysToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
      % ---------------------- FOR EACH DAY
        
        
        for j=DaysToPlot;
            [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
            title(['day ' num2str(j)]);
            
            Ymean_all=[];
            Ymean_AFP_all=[];
            
            for syl=SylsUnique;
                syl=syl{1};
                
                
                Ymean_AFP=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl).Pitch_PBSMinusMusc_MinusBase(j);
                Ymean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow(j);
                
                Ymean_all=[Ymean_all Ymean];
                Ymean_AFP_all=[  Ymean_AFP_all Ymean_AFP];
                    
            end
                        X=1:length(Ymean_all);
                lt_plot_bar(X-0.1, Ymean_all, {'Color','w', 'BarWidth',0.25});
                lt_plot_bar(X+0.1, Ymean_AFP_all, {'Color','b','BarWidth',0.25});
           
                set(gca, 'Xtick', X);
                set(gca, 'XTickLabel', SylsUnique);
                rotateXLabels(gca, 90);
        end
        
        lt_subtitle([birdname '-' exptname ', learn(wh), AFPbias(bl)']);
        
    end
end

            

%% BIN DAYS AND PLOT LEARNING AND AFP BIAS AT TARGET VS OTHER SYL TYPES

% methods to use to create day bins
divide_into_two=1;
% ExperimentsToExclude={'pu53wh88-SeqDepPitchLMAN','pu11wh87-SeqDepPitchLMAN3','gr41gr90-SeqDepPitchLMAN','gr41gr90-SeqDepPitchLMAN2'};
ExperimentsToExclude={'pu53wh88-SeqDepPitchLMAN'};
exclude_experiments =1;
clear OUTPUT;

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % +++++++++++++++++++++++ 0) Throw out if excluded
        if exclude_experiments==1;
            if any(strcmp([birdname '-' exptname], ExperimentsToExclude))
                continue;
            end
        end
        
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        % +++++++++++ 1) FIGURE OUT WHAT DAY BINS TO USE
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        
        if divide_into_two==1; % halfway between 1st and last musc day, not based on days of inactivaiton within.  inclusive for early.
            DayBinsToPlot={};
            musc_days_within_window=[];
            
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
                WNdayLast=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1)-1;
            else
                WNdayLast=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            musc_days_within_window=days_with_musc(days_with_musc>=WNday1 & days_with_musc<=WNdayLast);
            halfway_day=musc_days_within_window(1)+ceil(musc_days_within_window(end)-musc_days_within_window(1))/2;
            
            DayBinsToPlot{1}=musc_days_within_window(musc_days_within_window<=halfway_day);
            DayBinsToPlot{2}=musc_days_within_window(musc_days_within_window>halfway_day);
        elseif (0)
            
        end
        
        disp([birdname '-' exptname ' day bins: '  ]);
        cellfun(@disp, DayBinsToPlot);
        disp(' ');
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        % +++++++++++ 2) COLLECT DATA
        if sum(~cellfun(@isempty, DayBinsToPlot))<2;
            % then have only 0 or 1 bin with data
            continue;
        end
            
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;        
        for syl=SylsUnique;
            syl=syl{1};
            
            % ---------------------- FOR BIN
            for j=1:length(DayBinsToPlot);
                
                DaysToPlot=DayBinsToPlot{j};
                
                FFvals_learn_all=[];
                FFvals_musc_all=[];
                for daynum=DaysToPlot;
                    % ---- COLLECT RAW VALULES FOR THIS SYL (for this day)
                    
                    FFvals_Learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{daynum};
                    FFvals_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{daynum};
                    
                    % ------
                    FFvals_learn_all=[FFvals_learn_all FFvals_Learning];
                    FFvals_musc_all=[FFvals_musc_all FFvals_MUSC];
                end
                
                % ===== OUTPUT: for this daybin, add data for one more syl
                % ------------- INITIATE
                if ~exist('OUTPUT','var');
                    x=length(DayBinsToPlot);
                OUTPUT(length(DayBinsToPlot)).Learning_mean=[];
                OUTPUT(length(DayBinsToPlot)).Learning_sem=[];
                
                OUTPUT(length(DayBinsToPlot)).MP_mean=[];
                OUTPUT(length(DayBinsToPlot)).MP_sem=[];
               
                OUTPUT(x).AFP_mean=[];
                
                OUTPUT(length(DayBinsToPlot)).Similar=[];
                OUTPUT(length(DayBinsToPlot)).Presimilar=[];
                OUTPUT(length(DayBinsToPlot)).Target=[];      
                OUTPUT(x).birdnum=[];
                OUTPUT(x).exptnum=[];
                OUTPUT(x).targ_learn_dir=[];
                end
                
                    
                % ------------------ COLLECT
                OUTPUT(j).Learning_mean=[OUTPUT(j).Learning_mean mean(FFvals_learn_all)];
                OUTPUT(j).Learning_sem=[OUTPUT(j).Learning_sem lt_sem(FFvals_learn_all)];
                
                OUTPUT(j).MP_mean=[OUTPUT(j).MP_mean mean(FFvals_musc_all)];
                OUTPUT(j).MP_sem=[OUTPUT(j).MP_sem lt_sem(FFvals_musc_all)];
               
                OUTPUT(j).AFP_mean=[OUTPUT(j).AFP_mean mean(FFvals_learn_all)-mean(FFvals_musc_all)];
                
                OUTPUT(j).Similar=[OUTPUT(j).Similar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                OUTPUT(j).Presimilar=[OUTPUT(j).Presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
                OUTPUT(j).Target=[OUTPUT(j).Target SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                OUTPUT(j).birdnum=[OUTPUT(j).birdnum i];
                OUTPUT(j).exptnum=[OUTPUT(j).exptnum ii];
                
                OUTPUT(j).targ_learn_dir=[OUTPUT(j).targ_learn_dir SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir];
                
            end
        end
    end
end



%% -- backup of above
% methods to use to create day bins
divide_into_two=1;
ExperimentsToExclude={'pu53wh88-SeqDepPitchLMAN','pu11wh87-SeqDepPitchLMAN3','gr41gr90-SeqDepPitchLMAN','gr41gr90-SeqDepPitchLMAN2'};
exclude_experiments =1;
clear OUTPUT;

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % +++++++++++++++++++++++ 0) Throw out if excluded
        if exclude_experiments==1;
            if any(strcmp([birdname '-' exptname], ExperimentsToExclude))
                continue;
            end
        end
        
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        % +++++++++++ 1) FIGURE OUT WHAT DAY BINS TO USE
        days_with_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        
        if divide_into_two==1; % halfway between 1st and last musc day, not based on days of inactivaiton within.  inclusive for early.
            DayBinsToPlot={};
            musc_days_within_window=[];
            
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
                WNdayLast=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds(1)-1;
            else
                WNdayLast=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            end
            musc_days_within_window=days_with_musc(days_with_musc>=WNday1 & days_with_musc<=WNdayLast);
            halfway_day=musc_days_within_window(1)+ceil(musc_days_within_window(end)-musc_days_within_window(1))/2;
            
            DayBinsToPlot{1}=musc_days_within_window(musc_days_within_window<=halfway_day);
            DayBinsToPlot{2}=musc_days_within_window(musc_days_within_window>halfway_day);
        elseif (0)
            
        end
        
        disp([birdname '-' exptname ' day bins: '  ]);
        cellfun(@disp, DayBinsToPlot);
        disp(' ');
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        % +++++++++++ 2) COLLECT DATA
        if sum(~cellfun(@isempty, DayBinsToPlot))<2;
            % then have only 0 or 1 bin with data
            continue;
        end
            
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;        
        for syl=SylsUnique;
            syl=syl{1};
            
            % ---------------------- FOR BIN
            for j=1:length(DayBinsToPlot);
                
                DaysToPlot=DayBinsToPlot{j};
                
                FFvals_learn_all=[];
                FFvals_musc_all=[];
                for daynum=DaysToPlot;
                    % ---- COLLECT RAW VALULES FOR THIS SYL (for this day)
                    
                    FFvals_Learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{daynum};
                    FFvals_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{daynum};
                    
                    % ------
                    FFvals_learn_all=[FFvals_learn_all FFvals_Learning];
                    FFvals_musc_all=[FFvals_musc_all FFvals_MUSC];
                end
                
                % ===== OUTPUT: for this daybin, add data for one more syl
                % ------------- INITIATE
                if ~exist('OUTPUT','var');
                    x=length(DayBinsToPlot);
                OUTPUT(length(DayBinsToPlot)).Learning_mean=[];
                OUTPUT(length(DayBinsToPlot)).Learning_sem=[];
                
                OUTPUT(length(DayBinsToPlot)).MP_mean=[];
                OUTPUT(length(DayBinsToPlot)).MP_sem=[];
               
                OUTPUT(x).AFP_mean=[];
                
                OUTPUT(length(DayBinsToPlot)).Similar=[];
                OUTPUT(length(DayBinsToPlot)).Presimilar=[];
                OUTPUT(length(DayBinsToPlot)).Target=[];      
                OUTPUT(x).birdnum=[];
                OUTPUT(x).exptnum=[];
                OUTPUT(x).targ_learn_dir=[];
                end
                
                    
                % ------------------ COLLECT
                OUTPUT(j).Learning_mean=[OUTPUT(j).Learning_mean mean(FFvals_learn_all)];
                OUTPUT(j).Learning_sem=[OUTPUT(j).Learning_sem lt_sem(FFvals_learn_all)];
                
                OUTPUT(j).MP_mean=[OUTPUT(j).MP_mean mean(FFvals_musc_all)];
                OUTPUT(j).MP_sem=[OUTPUT(j).MP_sem lt_sem(FFvals_musc_all)];
               
                OUTPUT(j).AFP_mean=[OUTPUT(j).AFP_mean mean(FFvals_learn_all)-mean(FFvals_musc_all)];
                
                OUTPUT(j).Similar=[OUTPUT(j).Similar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                OUTPUT(j).Presimilar=[OUTPUT(j).Presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
                OUTPUT(j).Target=[OUTPUT(j).Target SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                OUTPUT(j).birdnum=[OUTPUT(j).birdnum i];
                OUTPUT(j).exptnum=[OUTPUT(j).exptnum ii];
                
                OUTPUT(j).targ_learn_dir=[OUTPUT(j).targ_learn_dir SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir];
                
            end
        end
    end
end


%% PLOT - GROSS AVERAGES
count=1;
SubplotsPerFig=12;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


% ============ Targets
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Pitch shift');
inds=OUTPUT(1).Target==1;
Ymean_all=[];
Ysem_all=[];
Yvals_all={};

% ----------------- LEARNING
Ymean=[];
Ysem=[];
Yvals={};

for j=1:length(OUTPUT);
    Yvals{j}=OUTPUT(j).Learning_mean(inds).*(OUTPUT(j).targ_learn_dir(inds)); % normalize to targ learn dir
    Ymean(j)=mean(Yvals{j});
    Ysem(j)=lt_sem(Yvals{j});
end

Ymean_all= [Ymean_all Ymean];
Ysem_all=[Ysem_all Ysem];
Yvals_all=[Yvals_all Yvals];


% ----------------------- AFP BIAS
Ymean=[];
Ysem=[];
Yvals={};

for j=1:length(OUTPUT);
    Yvals{j}=OUTPUT(j).AFP_mean(inds).*(OUTPUT(j).targ_learn_dir(inds)); % normalize to targ learn dir
    Ymean(j)=mean(Yvals{j});
    Ysem(j)=lt_sem(Yvals{j});
end

Ymean_all= [Ymean_all Ymean];
Ysem_all=[Ysem_all Ysem];
Yvals_all=[Yvals_all Yvals];


% -------------------------------- PLOT
% --- plot raw
for j=1:length(Yvals_all);
plot(j, Yvals_all{j}, 'ob', 'MarkerSize',5);
end
% --- plot means
Xlabels={'LEARNING-early','LEARNING-late','AFPbias-early','AFPbias-late'};
X=1:length(Ymean_all);
lt_plot_bar(X, Ymean_all, {'Errors',Ysem_all});

set(gca, 'XTick',X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45);




% ============================== AFP BIAS as fraction of learning
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('AFP bias as fraction of learning');
ylabel('fraction (paired divisions]');


Ymean=[];
Ysem=[];

for j=1:length(OUTPUT);
    Yvals_Learn=OUTPUT(j).Learning_mean(inds).*(OUTPUT(j).targ_learn_dir(inds)); % normalize to targ learn dir
    Yvals_AFP=OUTPUT(j).AFP_mean(inds).*(OUTPUT(j).targ_learn_dir(inds)); % normalize to targ learn dir
   
    % paired fractions
    Yvals_AFPfractionOfLearn=Yvals_AFP./Yvals_Learn;
    
    Ymean(j)=mean(Yvals_AFPfractionOfLearn);
    Ysem(j)=lt_sem(Yvals_AFPfractionOfLearn);
end



% -------------------------------- PLOT
% --- plot means
Xlabels={'early', 'late'};
X=1:length(Ymean);
lt_plot_bar(X, Ymean, {'Errors',Ysem});

set(gca, 'XTick',X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45);


                    
                    
                