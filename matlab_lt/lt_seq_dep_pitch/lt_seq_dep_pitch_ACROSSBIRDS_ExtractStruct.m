%% LT 7/21/15 - give it SeqDepPitch_AcrossBirds structure. Tell it how to filter. It gives you filtered struct

% e.g.
% filter = 'LMAN'   gives you only experiments with LMAN inactivation.
% filter = 'good_learning' only if learning is pass 100hz. (also throws out
%   if learning metric is nan
% filter = 'notLMAN' removes LMAN experiments
% filter = 'learnig_metric'; removed experiments where learning metric is
% nan;


function [SeqDepPitch_AcrossBirds_FILTERED, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter)

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% LEARNING METRIC -
if strcmp(filter, 'learning_metric');
    
        % copy strcuture, save backup.
    if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
        SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
    end
       
    
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
            end
        end
    end
    
end

%% LMAN - SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
if strcmp(filter, 'LMAN');
    
    % copy strcuture, save backup.
    if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
        SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
    end
    
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
    % copy strcuture, save backup.
    if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
        SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
    end
    
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
    
    % copy strcuture, save backup.
    if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
        SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
    end
    
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


