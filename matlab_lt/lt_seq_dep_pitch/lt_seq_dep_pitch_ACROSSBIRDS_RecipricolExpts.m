function [OUTSTRUCT] = lt_seq_dep_pitch_ACROSSBIRDS_RecipricolExpts(SeqDepPitch_AcrossBirds, DoDiffType)
%% lt 5/24/17 - number fo contexts greater, "drags" down learning in the target context?

NumBirds = length(SeqDepPitch_AcrossBirds.birds);
minLearnThresh = 50; % in hz, to keep expt.
% DoDiffType = 1; % if 1, then compares to diff type. if 0, then to same type.

%% === for each experiment, extract information to plot

OUTSTRUCT= struct;
count = 1;
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % === FOR EACH SAME-TYPE SYL, SEARCH TO SEE IF THERE IS AN EXPT
        % WHERE THIS SYL WAS THE TARG SYL
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        SameTypeSyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];

        DiffTypeSyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS];
        
        if DoDiffType==1;
            TypeToCompare = DiffTypeSyls;
        else
        TypeToCompare = SameTypeSyls;
        end
        
        for j=1:length(TypeToCompare)
            syl = TypeToCompare{j};
            
            % -- go thru all experiments, find one that has this syl as
            % target
            for jj = (ii+1):numexperiments
                if ~strcmp(syl, SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.INFORMATION.targsyl)
                    continue
                end
                
                % make sure that the targsyl from previous expt exists as a
                % syl in this expt
                if ~any(strcmp(targsyl, SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.INFORMATION.SylFields_Unique));
                    continue
                end
                
                % --- GOOD! save information about this pair
                generalizationExpt1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
                generalizationExpt2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean_rel_targ;
                targsylExpt1 = targsyl;
                targsylExpt2 = syl;
                if datenum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay, 'ddmmmyyyy') < ...
                        datenum(SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.Data_PlotLearning.Params.SeqFilter.FirstDay, 'ddmmmyyyy')
                EarlierExpt = 1;
                else
                    EarlierExpt = 2;
                end
                
                % ---- ONLY SAVE IF BOTH EXPERIMENTS HAD DECENT LEARNING
                learningExpt1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
                learningExpt1 = learningExpt1 * SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                learningExpt2 = SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                learningExpt2 = learningExpt2 * SeqDepPitch_AcrossBirds.birds{i}.experiment{jj}.INFORMATION.targ_learn_dir;
                
                if EarlierExpt == 2; % then flip around
                OUTSTRUCT(count).generalizationExpt2 = generalizationExpt1;
                OUTSTRUCT(count).generalizationExpt1 = generalizationExpt2;
                OUTSTRUCT(count).targsylExpt2 = targsylExpt1;
                OUTSTRUCT(count).targsylExpt1 = targsylExpt2;
                OUTSTRUCT(count).EarlierExpt = 1;
                OUTSTRUCT(count).LearningExpt1and2 = [learningExpt2; learningExpt1];
                else
                 OUTSTRUCT(count).generalizationExpt1 = generalizationExpt1;
                OUTSTRUCT(count).generalizationExpt2 = generalizationExpt2;
                OUTSTRUCT(count).targsylExpt1 = targsylExpt1;
                OUTSTRUCT(count).targsylExpt2 = targsylExpt2;
                OUTSTRUCT(count).EarlierExpt = EarlierExpt;
                OUTSTRUCT(count).LearningExpt1and2 = [learningExpt1; learningExpt2];                   
                end
                count = count+1;
                
                % temp
%                 if generalizationExpt1<-0.3
%                     keyboard
%                 end
                
            end
            
        end
        
    end
end

%% 
% =========== 1) generalization
lt_figure; hold on;
xlabel('generalization (earlier expt)');
ylabel('generalization (later expt)');
title('recipricol expts (e.g. A-B, then B-A)');
% -- earlier expt is expt 1
assert(all([OUTSTRUCT.EarlierExpt] ==1), 'PROBLEM - was assuming that all earlier expt occur first in SeqDepPitch');


% -- only keep those pass minimum laerning thresh
inds = min([OUTSTRUCT.LearningExpt1and2]) > minLearnThresh;

X = [OUTSTRUCT.generalizationExpt1];
Y = [OUTSTRUCT.generalizationExpt2];
lt_plot_45degScatter(X(inds), Y(inds), 'r', 1);
lt_regress(Y(inds), X(inds) , 1, 0, 1, 1, 'r', 0);

xlim([-0.6 1.4]); ylim([-0.6 1.4]);

% =========== 1) learning
lt_figure; hold on;
xlabel('learning (hz) (earlier expt)');
ylabel('learning (hz) (later expt)');

% -- earlier expt is expt 1
assert(all([OUTSTRUCT.EarlierExpt] ==1), 'PROBLEM - was assuming that all earlier expt occur first in SeqDepPitch');

tmp = [OUTSTRUCT.LearningExpt1and2];
X = tmp(1, :);
Y = tmp(2,:);
lt_plot_45degScatter(X, Y, 'r', 1);
lt_regress(Y, X , 1, 0, 1, 1, 'r', 0);
xlim([-50 250]); ylim([-50 250]);


