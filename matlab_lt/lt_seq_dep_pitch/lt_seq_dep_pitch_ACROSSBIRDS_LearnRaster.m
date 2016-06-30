function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LearnRaster(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 1/11/16 -


%%

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% ALL EXPT, LEARNING, IN TARG ORDER
lt_figure; hold on;
xlabel('shift from baseline');

% ==== determine what y-value to plot each experiment at (based on learning
% magnitude at target)
LearningAll=[];
Ymax=0; % to tally num expts
for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexpts;

        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % === collect
        LearningAll=[LearningAll learning];
        Ymax=Ymax+1;
    end
end

% ----- SORT LEARNING in order
[~, inds]=sort(LearningAll);
Yvals=1:Ymax;
Yvals=Yvals(inds);



ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir; % flip if needed
            
            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
                drift=drift*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                % --- STSS
                if similar & presim;
                    lt_plot(learning, Y, {'Color','b'});
                elseif similar & ~presim
                    lt_plot(learning, Y, {'Color','c'});
                elseif ~similar & presim
                    lt_plot(learning, Y, {'Color','r'});
                elseif ~similar & ~presim
                    plot(learning, Y, 'om');
                end
%                                 if similar & presim;
%                     plot(learning, Y, 'ob');
%                 elseif similar & ~presim
%                     plot(learning, Y, 'oc');
%                 elseif ~similar & presim
%                     plot(learning, Y, 'or');
%                 elseif ~similar & ~presim
%                     plot(learning, Y, 'om');
%                 end

            end
               
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','k');
 

% ==== PUT LINES FOR BASELINE DRIFT
disp(['Num syls with drift: ' num2str(NumSylsWithDrift)]);
disp(['Num syls no drift: ' num2str(NumSylsNoDrift)]);

% -- 2.5 and 97.5%ile [ALL]
DriftBounds=prctile(DriftBaseAll, [2.5 97.5]);
line([DriftBounds(1) DriftBounds(1)], ylim, 'Color','k','LineStyle','--');
line([DriftBounds(2) DriftBounds(2)], ylim, 'Color','k','LineStyle','--');
disp(['[ALL] 2.5 and 97.5% tile = ' num2str(DriftBounds(1)) ' & ' num2str(DriftBounds(2))]);

lt_plot_annotation(1, 'dashed line: 2.5/97.5 %-ile of drift of subset of all syls', 'k');

% % -- 2.5 and 97.5%ile [TARG]
% DriftBoundsTarg=prctile(DriftBaseTarg, [2.5 97.5]);
% line([DriftBoundsTarg(1) DriftBoundsTarg(1)], ylim);
% line([DriftBoundsTarg(2) DriftBoundsTarg(2)], ylim);
% disp(['[TARG] 2.5 and 97.5% tile = ' num2str(DriftBoundsTarg(1)) ' & ' num2str(DriftBoundsTarg(2))]);
% 

%% ==== output drift bounds


PARAMS.baselineDrift.AllDriftVals=DriftBaseAll;
PARAMS.baselineDrift.AllDriftVals_zscore=DriftBaseAll;
PARAMS.baselineDrift.NumSylsWithDrift=NumSylsWithDrift;
PARAMS.baselineDrift.NumSylsNoDrift=NumSylsNoDrift;
PARAMS.baselineDrift.prctiles_2_5and97_5=[DriftBounds(1) DriftBounds(2)];



%% PLOT AS ABOVE BUT DO NOT COLOR NONTARGETS

% ==== determine what y-value to plot each experiment at (based on learning
% magnitude at target)
LearningAll=[];
Ymax=0; % to tally num expts
for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexpts;

        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % === collect
        LearningAll=[LearningAll learning];
        Ymax=Ymax+1;
    end
end

% ----- SORT LEARNING in order
[~, inds]=sort(LearningAll);
Yvals=1:Ymax;
Yvals=Yvals(inds);


% ============ PLOT CDF OF LEARNING
lt_figure; hold on;
title('cdf of learning at targ');
lt_plot_cdf(LearningAll);


lt_figure; hold on;
xlabel('shift from baseline');

ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir; % flip if needed
            
            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
                drift=drift*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                    plot(learning, Y, 'ok');
            end
               
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','b');
 

% ==== PUT LINES FOR BASELINE DRIFT
disp(['Num syls with drift: ' num2str(NumSylsWithDrift)]);
disp(['Num syls no drift: ' num2str(NumSylsNoDrift)]);

% -- 2.5 and 97.5%ile [ALL]
DriftBounds=prctile(DriftBaseAll, [2.5 97.5]);
line([DriftBounds(1) DriftBounds(1)], ylim, 'Color','b','LineStyle','--');
line([DriftBounds(2) DriftBounds(2)], ylim, 'Color','b','LineStyle','--');
disp(['[ALL] 2.5 and 97.5% tile = ' num2str(DriftBounds(1)) ' & ' num2str(DriftBounds(2))]);

lt_plot_annotation(1, 'dashed line: 2.5/97.5 %-ile of baseline drift of subset of all syls', 'k');


%% PLOT AS ABOVE but just color same type and diff type
lt_figure; hold on;
xlabel('shift from baseline');

% ==== determine what y-value to plot each experiment at (based on learning
% magnitude at target)
LearningAll=[];
Ymax=0; % to tally num expts
for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexpts;

        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % === collect
        LearningAll=[LearningAll learning];
        Ymax=Ymax+1;
    end
end

% ----- SORT LEARNING in order
[~, inds]=sort(LearningAll);
Yvals=1:Ymax;
Yvals=Yvals(inds);



ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning=learning*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir; % flip if needed
            
            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_mean;
                drift=drift*SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                % --- STSS
                if similar
                    lt_plot(learning, Y, {'Color','b'});
                elseif ~similar
                    lt_plot(learning, Y, {'Color','r'});
                end
%                                 if similar & presim;
%                     plot(learning, Y, 'ob');
%                 elseif similar & ~presim
%                     plot(learning, Y, 'oc');
%                 elseif ~similar & presim
%                     plot(learning, Y, 'or');
%                 elseif ~similar & ~presim
%                     plot(learning, Y, 'om');
%                 end

            end
               
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','k');
 

% ==== PUT LINES FOR BASELINE DRIFT
disp(['Num syls with drift: ' num2str(NumSylsWithDrift)]);
disp(['Num syls no drift: ' num2str(NumSylsNoDrift)]);

% -- 2.5 and 97.5%ile [ALL]
DriftBounds=prctile(DriftBaseAll, [2.5 97.5]);
line([DriftBounds(1) DriftBounds(1)], ylim, 'Color','k','LineStyle','--');
line([DriftBounds(2) DriftBounds(2)], ylim, 'Color','k','LineStyle','--');
disp(['[ALL] 2.5 and 97.5% tile = ' num2str(DriftBounds(1)) ' & ' num2str(DriftBounds(2))]);

lt_plot_annotation(1, 'dashed line: 2.5/97.5 %-ile of drift of subset of all syls', 'k');

% % -- 2.5 and 97.5%ile [TARG]
% DriftBoundsTarg=prctile(DriftBaseTarg, [2.5 97.5]);
% line([DriftBoundsTarg(1) DriftBoundsTarg(1)], ylim);
% line([DriftBoundsTarg(2) DriftBoundsTarg(2)], ylim);
% disp(['[TARG] 2.5 and 97.5% tile = ' num2str(DriftBoundsTarg(1)) ' & ' num2str(DriftBoundsTarg(2))]);
% 


%% ==== PLOT ORDERED AS ABOVE, BUT PLOT GENERALIZATION, NOT RAW SHIFT [order of absolute learning magn]

lt_figure; hold on;
title('generalization');
xlabel('generalization');


ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;

            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                % --- STSS
                if similar & presim;
                    lt_plot(learning, Y, {'Color','b'});
                elseif similar & ~presim
                    lt_plot(learning, Y, {'Color','c'});
                elseif ~similar & presim
                    lt_plot(learning, Y, {'Color','r'});
                elseif ~similar & ~presim
                    plot(learning, Y, 'om');
                end
%                                 if similar & presim;
%                     plot(learning, Y, 'ob');
%                 elseif similar & ~presim
%                     plot(learning, Y, 'oc');
%                 elseif ~similar & presim
%                     plot(learning, Y, 'or');
%                 elseif ~similar & ~presim
%                     plot(learning, Y, 'om');
%                 end

            end
               
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','k');


% % ==== PUT LINES FOR BASELINE DRIFT
% disp(['Num syls with drift: ' num2str(NumSylsWithDrift)]);
% disp(['Num syls no drift: ' num2str(NumSylsNoDrift)]);
% 
% % -- 2.5 and 97.5%ile [ALL]
% DriftBounds=prctile(DriftBaseAll, [2.5 97.5]);
% line([DriftBounds(1) DriftBounds(1)], ylim, 'Color','k','LineStyle','--');
% line([DriftBounds(2) DriftBounds(2)], ylim, 'Color','k','LineStyle','--');
% disp(['[ALL] 2.5 and 97.5% tile = ' num2str(DriftBounds(1)) ' & ' num2str(DriftBounds(2))]);
% 
% lt_plot_annotation(1, 'dashed line: 2.5/97.5 %-ile of drift of subset of all syls', 'k');
% 
%     

%% ==== PLOT ORDERED AS ABOVE, BUT PLOT GENERALIZATION, NOT RAW SHIFT [order of absolute learning magn]
% [just same type, diff type]

lt_figure; hold on;
title('generalization');
xlabel('generalization');


ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;

            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                % --- STSS
                if similar
                    lt_plot(learning, Y, {'Color','b'});
                elseif ~similar
                    lt_plot(learning, Y, {'Color','r'});
                end
%                                 if similar & presim;
%                     plot(learning, Y, 'ob');
%                 elseif similar & ~presim
%                     plot(learning, Y, 'oc');
%                 elseif ~similar & presim
%                     plot(learning, Y, 'or');
%                 elseif ~similar & ~presim
%                     plot(learning, Y, 'om');
%                 end

            end
               
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','k');


% % ==== PUT LINES FOR BASELINE DRIFT
% disp(['Num syls with drift: ' num2str(NumSylsWithDrift)]);
% disp(['Num syls no drift: ' num2str(NumSylsNoDrift)]);
% 
% % -- 2.5 and 97.5%ile [ALL]
% DriftBounds=prctile(DriftBaseAll, [2.5 97.5]);
% line([DriftBounds(1) DriftBounds(1)], ylim, 'Color','k','LineStyle','--');
% line([DriftBounds(2) DriftBounds(2)], ylim, 'Color','k','LineStyle','--');
% disp(['[ALL] 2.5 and 97.5% tile = ' num2str(DriftBounds(1)) ' & ' num2str(DriftBounds(2))]);
% 
% lt_plot_annotation(1, 'dashed line: 2.5/97.5 %-ile of drift of subset of all syls', 'k');
% 
%     
  
%% PLOT GENERALIZATION, BUT DON'T COLOR NONTARGS
lt_figure; hold on;
title('generalization');
xlabel('generalization');


ExptCount=1;
% ====== COLLECT ALL EXPERIMENTS
BirdnamesAll={};
ExptnamesAll={};
YtickLabAll={};
DriftBaseAll=[];
NumSylsWithDrift=0;
NumSylsNoDrift=0;
DriftBaseTarg=[];

for i=1:NumBirds
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        Y=find(Yvals==ExptCount);
        
            BirdnamesAll{Y}=birdname;
            ExptnamesAll{Y}=exptname;
            YtickLabAll{Y}=[birdname(1:4) '-' exptname(end-4:end)];
            
            
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
                
            
            % --- collect learning
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;

            % --- status
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % === BASELINE DRIFT COLLECT
            if isfield (SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'DRIFT_UsingTwoDays');
                drift=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).DRIFT_UsingTwoDays.zscore_rel_targ;
                
                DriftBaseAll=[DriftBaseAll drift];
                
                % --- only targets
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    DriftBaseTarg=[DriftBaseTarg drift];
                end
                
                NumSylsWithDrift=NumSylsWithDrift+1;
            else
                NumSylsNoDrift=NumSylsNoDrift+1;
            end
            
            % ==== plot
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                lt_plot(learning, Y);
            else
                    plot(learning, Y, 'ok');
            end
            
            
        end
        ExptCount=ExptCount+1;
    end
end

            
          
set(gca, 'Ytick', 1:ExptCount-1);
set(gca, 'YtickLabel', YtickLabAll);

line([0 0], ylim, 'Color','k');

