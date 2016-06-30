%% LT 7/21/15 - give it SeqDepPitch_AcrossBirds structure. Tell it how to filter. It gives you filtered struct

% e.g.
% filter = 'LMAN'   gives you only experiments with LMAN inactivation.
% filter = 'good_learning' only if learning is pass 100hz. (also throws out
%   if learning metric is nan
% filter = 'notLMAN' removes LMAN experiments
% filter = 'learnig_metric'; removed experiments where learning metric is
% nan;
% filter = 'ThrowOutIfTargInPosTwoOrHigherInRepeat' 
% filter = 'ThrowOutIfTargInPosThreeOrHigherInRepeat'
% filter = 'RepeatExptsOnly';
% filter = 'MultiDirAreInRepeats'
% filter = 'RemoveStartTwoTargs'
% filter = 'AdHocExptRemove'


function [SeqDepPitch_AcrossBirds_FILTERED, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter)

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% filter = DelayExptStart
% i.e. did not sing on first 1+ day once WN start

if strcmp(filter, 'DelayExptStart');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode>0
                % then had delay - throw out.
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
end


%% filter = learning_range

if strcmp(filter, 'learning_range');
    %     threshold=0.8; % zscore
    threshold_min=input('what is min zscore threshold for learning? ');
    threshold_max=input('what is max zscore thresh for learning? ');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            if  isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
                % - no learning value
                expts_to_remove{i}=[expts_to_remove{i} ii];
            else
                learnmetr=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
                
                % -- flip sign if negative
                targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                learnmetr=learnmetr*targdir;
                
                if learnmetr<threshold_min | learnmetr>threshold_max
                    % then throw out
                    expts_to_remove{i}=[expts_to_remove{i} ii];
                end
                
            end
        end
    end
    
end

%% filter = 'AdHocExptRemove'
% remove expt which have poor learning even late, or other abnormalities.
% diff from learning filter because later is based on day 3-4.
if strcmp(filter, 'AdHocExptRemove')
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            if strcmp(birdname, 'pu35wh17') & strcmp(exptname, 'SeqDepPitch3');
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
            if strcmp(birdname, 'rd28pu64') & strcmp(exptname, 'SeqDepPitchLMAN2');
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
            
        end
    end
end


%% filter = 'RemoveStartTwoTargs'

if strcmp(filter, 'RemoveStartTwoTargs')
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            numtargs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs;
            
            if numtargs==2;
                % mult targets. note down index and remove later
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
        
    end
    
    
    
end


%% MultiDirAreInRepeats
if strcmp(filter, 'MultiDirAreInRepeats');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            % === check whether there is repeats data
            if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}, 'DataRepeats');
                % then throw out
                expts_to_remove{i}=[expts_to_remove{i} ii];
            elseif ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION, 'MultiDirSyls');
                % if no multidir, then throw out
                expts_to_remove{i}=[expts_to_remove{i} ii];
            else
                % if there are repeats and multidir, make sure the
                % multidir are in repeats
                MultidirSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
                
                for j=1:length(MultidirSyls);
                    syl=MultidirSyls{j};
                    
                    if ~any(strcmp(syl, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.SylsInRepeat));
                        % then this syl is not part of targ syl repeats
                            expts_to_remove{i}=[expts_to_remove{i} ii];
                        break
                    end
                end
                
            end
            
            
        end
        
    end
end



%% RepeatExptsOnly

if strcmp(filter, 'RepeatExptsOnly');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            % === check whether there is repeats data
            if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}, 'DataRepeats');
                % then throw out
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
        
    end
end


%% Throw out if target is in repeat (pos 2+)

if strcmp(filter, 'ThrowOutIfTargInPosTwoOrHigherInRepeat');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            % if presyl is the same syl type as the target,
            % then this experiment must be targeting repeats
            presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).preceding_syl;
            
            if isnan(presyl)
                % then cannot be in pos 2+ of repeats
                continue
            end
            
            % CONVERT TO SINGLE SYL (lower)
            if length(presyl)>1
                [presyl]=regexp(presyl,'[A-Z]', 'match');
                presyl=lower(presyl);
            else
                presyl=lower(presyl);
            end
            
            % COMPARE PRESYL TO TARGET
            target_singlesyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;
            
            if strcmp(presyl, target_singlesyl);
                % then throw out
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
            %                % sanity check
            %             if regexp(presyl, '[A-Z]');
            %                 disp(presyl);
            %             end
            %
            %             if length(presyl)>1
            %                 disp(presyl)
            %             end
            %
        end
        
    end
    
end

%% Throw out if target is in repeat (pos 3+)

if strcmp(filter, 'ThrowOutIfTargInPosThreeOrHigherInRepeat');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            % collect 2 syl back and one syl back
            two_syl_back=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).two_syl_back;
            
            if isnan(two_syl_back);
                % then can't be pos 2 in repeat,
                continue
            end
            
            two_syl_back=lower(two_syl_back);
            assert(length(two_syl_back)==1, 'problem - two syl back is more than one syl');

            presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).preceding_syl;            
            if isnan(presyl)
                % then cannot be in pos 2+ of repeats
                continue
            end

            presyl=lower(presyl);
            assert(length(presyl)==1, 'problem - presyl is more than one syl');

            
%             % CONVERT TO SINGLE SYL (lower)
%             if length(presyl)>1
%                 [presyl]=regexp(presyl,'[A-Z]', 'match');
%                 presyl=lower(presyl);
%             else
%                 presyl=lower(presyl);
%             end
%             
%             
            % GET TARGET SINGLE SYL
            target_singlesyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;
            
            
            % ==== compare target to presyl to 2sylback
            if strcmp(presyl, target_singlesyl) & strcmp(presyl, two_syl_back)
                % then throw out
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
            
        end
        
    end
    
end


%% DID NOT SING 1- 1 day at start
if strcmp(filter, 'did_not_sing');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            if  isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals{WNonInd})
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
end

%% DID NOT SING 2- two days at start
if strcmp(filter, 'did_not_sing2');
    
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            
            if  isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals{WNonInd}) ...
                    & isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals{WNonInd+1})
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
end


%% LEARNING METRIC -
if strcmp(filter, 'learning_metric');
%     threshold=0.8; % zscore
    threshold=input('what is zscore threshold for learning? ');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            if  isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
                expts_to_remove{i}=[expts_to_remove{i} ii];
            elseif abs(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean)<threshold;
                expts_to_remove{i}=[expts_to_remove{i} ii];
                
            end
        end
    end
    
end

%% LMAN - SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
if strcmp(filter, 'LMAN');
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0;
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
    
end

%% GOOD LEARNING - 'good_learning' only if learning is pass 100hz.
if strcmp(filter, 'good_learning');
    learning_threshold=80;
    disp(['removing expts with < ' num2str(learning_threshold) 'hz learning']);
    
    %     % copy strcuture, save backup.
    %     if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
    %         SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
    %     end
    %
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            % ===== WHAT was learning at target?
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.consolid_start_meanFF_minbase;
            
            if abs(learning)<learning_threshold;
                expts_to_remove{i}=[expts_to_remove{i} ii];
            elseif  isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean);
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
    
end


%% NOT LMAN

if strcmp(filter, 'notLMAN');
    
    
    % remove experiments - first figure out what inds to remove
    expts_to_remove=[];
    for i=1:NumBirds;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        expts_to_remove{i}=[];
        for ii=1:NumExperiments;
            
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
                expts_to_remove{i}=[expts_to_remove{i} ii];
            end
        end
    end
    
    
end

%% END - REMOVE EXPTS, AND GENERAL THINGS

% REMOVE
disp(' ');
disp(['REMOVED EXPERIMENTS, WITH FILTER: ' filter]);


% actually remove them
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


% If any birds have nothing, then remove the birds
birdstoremove=[];
for i=1:NumBirds;
    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment)
        birdstoremove=[birdstoremove i];
    end
end

SeqDepPitch_AcrossBirds.birds(birdstoremove)=[];



% -----
SeqDepPitch_AcrossBirds_FILTERED=SeqDepPitch_AcrossBirds;
NumBirds=length(SeqDepPitch_AcrossBirds.birds);



% ---- display all birds remaining:
disp(' ');
disp('REMAINING BIRDS-EXPERIMENTS:');
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        disp([birdname ' - ' exptname]);
        
    end
end


