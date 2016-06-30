function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_EarlyDays(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 10/4/15 - some experiments, i.e. only for pu53 currently, show hard to understand early effects - large changes in pitch on fisrt few days, in opposite direction of learning.
% Throw those days out, so that z-score learning analysis automaticalyl
% starts counting number of days starting when actual data starts.

% Criteria for throwing out - shift by target is in opposite direction from
% learning, greater than a threshold.


%% Params
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

% days to check, WN1 to WN5 (max)
days_to_check_max=5; % WN day 5



%% Find experiemnts that pass criterion and remove those days
disp(' ');
disp('REMOVING BAD DAYS FOR: ');

for i=1:NumBirds;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % =============================== REMOVE STRUCTURE FIELDS THAT I DON'T UYSE THAT
        % MIGHT BE CONTRADICTIROY
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning, 'DataMatrix_Targ');
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning=rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning, 'DataMatrix_Targ');
        end
        
        % ================================= SEARCH TO SEE IF NEED TO REMOVE AERLY DAYS
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        dayInds_tocheck=WNday1:WNday1+days_to_check_max-1;
        
        
        
        % ----------------- 1) Start from WN day 1 --> day with targ ff in
        % opposite dir?  if yes, remove data, go to next day. if no, then
        % continue to next experiment
        learning_zscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_zscore;
        bad_day_inds=[];
        for j=1:length(dayInds_tocheck);
            dayind=dayInds_tocheck(j);
            
            learning=learning_zscore(dayind);
            
            % Go to next day if day has no data
            if isnan(learning);
                continue
            end
            
            % break for loop (over days), if day does not fail learning
            % criterion. if day fails, record the day and keep going in
            % loop.
            if targ_learn_dir==1;
                if learning<-PARAMS.EarlyDays.threshold;
                    % bad day, throw out
                    bad_day_inds=[bad_day_inds dayind];
                else
                    break
                end
            elseif targ_learn_dir==-1;
                if learning>PARAMS.EarlyDays.threshold;
                    % bad day, throw out
                    bad_day_inds=[bad_day_inds dayind];
                else
                    break
                end
            end
        end
        
        if isempty(bad_day_inds);
            continue;
        end
        
        % ======== FROM HERE, ONLY IF ACTUALLY FOUND BAD DAYS
        % =============================================================================
        % display to user the days found
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        
        disp([birdname '-' exptname ', bad days: ' num2str(bad_day_inds)]);
        
        % ---------------------------------------------------  What are the day Inds I should keep?
        dayInds_tokeep=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals);
        dayInds_tokeep(bad_day_inds)=[];
        
        
        % ==============================
        % Remove bad day inds
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            % ------------------------------------------ 1) Remove fields that I don't use, they
            % hinder the subsampling function
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl)= ...
                rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl), 'GeneralizationFrom');
            if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).CatchSongStatus)<max(dayInds_tokeep);
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl)=...
                    rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl), 'CatchSongStatus');
            end
            
            
            % --- Learning
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl)=...
                lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl), dayInds_tokeep);
            
            
            % -------- Hit rate
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays=...
                lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays, dayInds_tokeep);
            
            
            % -------- Remove regexpr data
            try % sometimes regexpr does not have enough days
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data(bad_day_inds)=[];
            catch err
            end
            
            
            
            
            
            % ------------------------------------------------------------ 3) Do same for LMAN data if it exists
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
                
                % -- remove things
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl)= ...
                    rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl), 'GeneralizationFrom');
                if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).CatchSongStatus)<max(dayInds_tokeep);
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl)=...
                        rmfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl), 'CatchSongStatus');
                end
                
                % --- Learning
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl)=...
                    lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl), dayInds_tokeep);
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl)=...
                    lt_structure_subsample_all_fields(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSCvsPBS.(syl), dayInds_tokeep);
                
                try % sometimes regexpr does not have enough days
                    
                    % ---- Remove regexpr data
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data_MUSC(bad_day_inds)=[];
                catch err
                end
                
            end
            
            
        end
    end
end


disp('done!');



