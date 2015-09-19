%% LT 4/13/15 - Compiles data and plots across experiments/birds.  Need to specific dir of saved structures
function lt_seq_dep_pitch_ACROSSBIRDS

%% TO DO:
% 1) plot similar and diff syls, separated by whether in same or diff motif


%% NEED DATA STRUCTURES
% from PlotLearning: AllDays_PlotLearning and Params
clear all; close all;

% Cell holding experiments {birdname, experient ID, NumTargets(e.g. 1, 2, 0(means 1 and more in diff epochs), directory of plotlearning, ConsolidationStart, ConsolidationEnd};
ExperimentList{1}={'pu53wh88','SyntaxDepPitchShift_abDowncbUp',2,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift/SeqFilterCompile',...
    '16Jul2014','21Jul2014'};
% ExperimentList{2}={'pu53wh88','SyntaxDepPitchShift_cbUp',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SyntaxDepPitchShift_cbUp/SeqFilterCompile',...
%     '04Aug2014','06Aug2014'}; THROW OUT AS MAX LEARNING ONLY 50HZ
ExperimentList{2}={'pu53wh88','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '06Nov2014','11Nov2014'};
ExperimentList{3}={'pu53wh88','SeqDepPitchShift2',1','/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '04Dec2014','11Dec2014'};
ExperimentList{4}={'pu53wh88','SeqDepPitchShift3',0,'/bluejay3/lucas/birds/pu53wh88/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '20Feb2015','24Feb2015'};

ExperimentList{5}={'pu11wh87','SyntaxDepPitchShift_abUP',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_abUP/SeqFilterCompile',...
    '12Jul2014','15Jul2014'};
ExperimentList{6}={'pu11wh87','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '03Dec2014','09Dec2014'};
ExperimentList{7}={'pu11wh87','SeqDepPitchShift3',0,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift3/SeqFilterCompile',...
    '26Feb2015','07Mar2015'};
ExperimentList{8}={'pu11wh87','SyntaxDepPitchShift_cbDOWN',1,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SyntaxDepPitchShift_cbDOWN/SeqFilterCompile',...
    '28Jul2014','03Aug2014'};
ExperimentList{9}={'pu11wh87','SeqDepPitchShift',2,'/bluejay3/lucas/birds/pu11wh87/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '02Nov2014','07Nov2014'};


ExperimentList{10}={'pu37wh20','SeqDepPitchShift',1,'/bluejay3/lucas/birds/pu37wh20/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '02Nov2014','03Nov2014'};

ExperimentList{11}={'pu64bk13','SeqDepPitchShift',0,'/bluejay3/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '17Dec2014','21Dec2014'};
ExperimentList{12}={'pu64bk13','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/pu64bk13/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '17Feb2015','22Feb2015'};

ExperimentList{13}={'gr41gr90','SeqDepPitchShift',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift/SeqFilterCompile',...
    '15Feb2015','24Feb2015'};
ExperimentList{14}={'gr41gr90','SeqDepPitchShift2',0,'/bluejay3/lucas/birds/gr41gr90/seq_dep_pitch_SeqDepPitchShift2/SeqFilterCompile',...
    '18Apr2015','23Apr2015'};

ExperimentList{15}={'rd23gr89','SeqDepPitch',1,'/bluejay4/lucas/birds/rd23gr89/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '09Jun2015','14Jun2015'};

ExperimentList{16}={'pu35wh17','SeqDepPitch',1,'/bluejay4/lucas/birds/pu35wh17/seq_dep_pitch_SeqDepPitch/SeqFilterCompile',...
    '10Jun2015','14Jun2015'};



%% EXTRACT EXPERIMENT DATA STRUCTURES
NumExpt=length(ExperimentList);
SeqDepPitch_AcrossBirds=struct('birds',[]); % holding structure

for i=1:length(ExperimentList);
    birdname=ExperimentList{i}{1};
    ExptID=ExperimentList{i}{2};
    NumTargs=ExperimentList{i}{3};
    PlotLearningDir=ExperimentList{i}{4};
    ConsolStartDate=ExperimentList{i}{5};
    ConsolEndDate=ExperimentList{i}{6};
    
    
    % Slot into structure
    if isempty(SeqDepPitch_AcrossBirds.birds);
        % then no birds slotted yet
        % start bird 1
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.PlotLearningDir=PlotLearningDir;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.NumTargs=NumTargs;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.ExptID=ExptID;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
        SeqDepPitch_AcrossBirds.birds{1}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
        SeqDepPitch_AcrossBirds.birds{1}.birdname=birdname;
        
        
    else
        % find the right bird to slot data into
        c=0;
        for j=1:length(SeqDepPitch_AcrossBirds.birds);
            if strcmp(SeqDepPitch_AcrossBirds.birds{j}.birdname,birdname);
                % then this is right bird
                % enter new data, new expt
                ind=length(SeqDepPitch_AcrossBirds.birds{j}.experiment)+1;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.PlotLearningDir=PlotLearningDir;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.NumTargs=NumTargs;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.ExptID=ExptID;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolStartDate=ConsolStartDate;
                SeqDepPitch_AcrossBirds.birds{j}.experiment{ind}.DATES.ConsolEndDate=ConsolEndDate;
                SeqDepPitch_AcrossBirds.birds{j}.birdname=birdname;
                % change c to show that bird found.
                c=1;
                
                continue % quit loop
            end
        end
        if c==0
            % then no bird was found, start new bird
            ind=length(SeqDepPitch_AcrossBirds.birds)+1;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.PlotLearningDir=PlotLearningDir;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.NumTargs=NumTargs;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.ExptID=ExptID;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolStartDate=ConsolStartDate;
            SeqDepPitch_AcrossBirds.birds{ind}.experiment{1}.DATES.ConsolEndDate=ConsolEndDate;
            SeqDepPitch_AcrossBirds.birds{ind}.birdname=birdname;
        end
    end
end



%-- LOAD AND COMPILE STRUCTURES
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        % load
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/Params']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params=tmp.Params;
        tmp=load([SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/AllDays_PlotLearning']);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning=tmp.AllDays_PlotLearning;
        
    end
end


% -- HOW many birds? 
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


% -- MAKE RANDOM EDITS (ad hoc)
 
% pu37wh20, % remove bccB as this was actually also a target (along with dccB);
% Don't need to, not in data nymore...
% for i=1:NumBirds;
%     if strcmp(SeqDepPitch_AcrossBirds.birds{i}.birdname,'pu37wh20');
%     
%         for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%             if strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
%                 
% %                 % remove stuff
% %                 ii=find(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsSame,'bccB'))
%                 
%             end
%     
%     end
% end


% CHECK DIRECTION OF PITCH SHIFT

%% PLOT SUMMARY OF METADATA OF EXPERIMENTS (E.G. SAMPLE SIZE)
% IN PROGRESS
% e..g bar plot, each col is a bird, with height as num experiments of
% different ytpes)




%% GET GLOBAL FEATURE PCA SPACE

% Summary of methods: Each bird contributes a syllable only a single time (even if
% performed multiple experiments). Take from just oen experiment
% (choose the middle one date-wise). each syllable contributes the same
% number of renditions (randomly choosing, the number size of smallest
% sample size). All randitions converted to a global z-score, and then
% PCA on that.


currdir=pwd;

% LOAD BASELINE FEATURE VECTORS FOR ALL RENDITIONS OF ALL UNIQUE SYLLABLES
% ACROSS ALL BIRDS

% ===1) LOAD StatsStruct if it has been made. If not, run program to make
% it.
for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    % ============= HOW MANY EXPERIMENTS?
    % if only one, then use this experiment's data
    % if more than one, then collect data from the middle experiment
    if numexperiments==1;
        expt_ind=1;
    else
        % take the middle experiment
        expt_ind=ceil(numexperiments/2);
    end
    
    
    % ============= WHAT ARE META DATA FROM THIS EXPERIMENT
    if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning, 'SylFields_Unique'); % unique syls
        syls_unique={}; % then collect unqiue syls from scratch
        for j=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder);
            syls_unique=[syls_unique SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{j}];
        end
        SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=syls_unique;
    end
    
    % Collect stuff
    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
    save_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.PlotLearningDir;
    
    
    % ============= COLLECT FEATURE VECTORS FOR ALL BASELINE RENDITIONS
    cd(save_dir)
    % --- Does saved structure already exist?
    % IS THE FOLDER STRUCTURED IN THE OLD OR NEW FORMAT?
    tmp=dir('AllDays_RawDatStruct.mat');
    if isempty(tmp)
        % Old format -->
        disp('PROBLEM - folder structured in old format - bad...');
    else
        % New format -->
        tmp=dir('DONE_StructureStats*');
        if ~isempty(tmp);
            % THEN SAVED STRUCTURE EXISTS - Load and extract data
            load AllDays_StructStatsStruct;
            load Params;
            
        else
            % THEN HAVE TO RUN CODE TO GET STATS STRUCT
            load Params;
            load AllDays_RawDatStruct;
            
            [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct);
            close all;
        end
    end
    cd(currdir);
    
    
    % ============ TAKE THAT FEATURE STRUCTURE AND SAVE IT TO THE GLOBAL
    % STRUCTURE
    SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.used_which_expt=expt_ind;
    SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.syls_unique=syls_unique;
    
    % extract data for each syl
    baseline_days=SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
    
    for k=1:length(syls_unique);
        syl=syls_unique{k};
        
        inds_from_baseline=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd<(baseline_days(end)+1);
        
        % ========= FINAL DATA
        featvect_tmp=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(inds_from_baseline,:);
        SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.legend=Params.StructureStats.FeatureLegend;
        SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data.(syl).fv=featvect_tmp;
    end
    
end
%

    
%% ====== PERFORM PCA ON ALL SYLLABLES
% == First collect all syllables into one large array (subsample to get same
% number of samples for each syllable)

% 1) ============= SAMPLE SIZE 
% First, display sample sizes:
disp('SAMPLE SIZES');
for i=1:NumBirds;
    disp(' ');
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    disp(['Bird '  num2str(i) ' (' birdname ')']);
    
    unique_syls=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data);
    for ii=1:length(unique_syls);
        syl=unique_syls{ii};
        
        N=size(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data.(syl).fv,1);
        
        disp([syl ': ' num2str(N)]);
    end
end

% Second, 
% Plot histogram of sample sizes:
N_all=[];
for i=1:NumBirds;
    unique_syls=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data);
    for ii=1:length(unique_syls);
        syl=unique_syls{ii};
        
        N=size(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data.(syl).fv,1);
        
        N_all=[N_all N];
    end
end

lt_figure; hold on;
title('histogram of sample sizes over all syls');
hist(N_all,50);

% Third, query user for sample size;
N_subsample=input('how many renditions to take for each syl? ');


% 2) ============== Subsample and collect all data
PCA.raw_data.FV_AllSylsAllBirds=[];

for i=1:NumBirds;

    unique_syls=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data);
    for ii=1:length(unique_syls);
        syl=unique_syls{ii};
    
        % SHUFFLE DATA and SUBSAMPLE
        N=size(SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data.(syl).fv,1);

        N_shuff=randperm(N);
        
        if N>=N_subsample;
            Inds=N_shuff(1:N_subsample);
        else
            Inds=N_shuff;
        end
        
        % --- Collect data
        featurevect_temp=SeqDepPitch_AcrossBirds.birds{i}.baseline_feature_vectors.data.(syl).fv(Inds,:);
        
        PCA.raw_data.FV_AllSylsAllBirds=[PCA.raw_data.FV_AllSylsAllBirds; featurevect_temp];
        
    end
end

PCA.raw_data.N_subsample=N_subsample;



% 3) =============== z-score all data
PCA.raw_data.FV_mean=mean(PCA.raw_data.FV_AllSylsAllBirds,1);
PCA.raw_data.FV_std=std(PCA.raw_data.FV_AllSylsAllBirds);

N=size(PCA.raw_data.FV_AllSylsAllBirds,1);

PCA.raw_data.FV_MinusMean=PCA.raw_data.FV_AllSylsAllBirds-repmat(PCA.raw_data.FV_mean,N,1);
PCA.raw_data.FV_zscore=PCA.raw_data.FV_MinusMean./repmat(PCA.raw_data.FV_std,N,1);



% 4) ================== PCA
[coeff, score, latent, tsquared, explained]=pca(PCA.raw_data.FV_zscore);

PCA.outputs.coeff=coeff; % orthonormal vector coefficients
PCA.outputs.score=score; % coords of original data in new coord system
PCA.outputs.explained=explained; % explained var
PCA.outputs.latent=latent; % variance explained of each PC
PCA.outputs.tsquared=tsquared; % hotelling t distance for each rend (i.e. var weighted distance from global centroid).

% ----------- PLOT THINGS ABOUT PCA OUTPUT
% Plot 1ST VS. 2ND PC, and 1st 3 in 3d
lt_figure; hold on;
plot3(score(:,1), score(:,2), score(:,3),'+');
title('PC3 vs. PC2 vs. PC1');

% DISPLAY PCA stats
lt_figure; hold on;
pareto(explained);
title('% variance explained by each PC');

% WHAT DO the PCs look like?
lt_figure; hold on;
title('PC coefficients (lines) for first 3 PCs, and all scores from data (normalized)');
biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',SeqDepPitch_AcrossBirds.birds{1}.baseline_feature_vectors.legend)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');



% 5) ================== FOR EACH SYL, get COM distance from target








%% SORT OUT ONLY THE EXPERIMENTS WITH ONLY ONE TARG, OR STARTED OUT WITH ONLY ONE TARG

% copy strcuture, the new one will have experiments removed that do no pass contingency
SeqDepPitch_AcrossBirds_ONETARG=SeqDepPitch_AcrossBirds;

% remove experiments - first figure out what inds to remove
expts_to_remove=[];
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment);
    
    expts_to_remove{i}=[];
    for ii=1:NumExperiments;
        numtargs=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.NumTargs;
        
        if numtargs>1;
            % mult targets. note down index and remove later
            expts_to_remove{i}=[expts_to_remove{i} ii];
        end
    end
end

% actually remove them
for i=1:length(expts_to_remove);
    if ~isempty(expts_to_remove{i});
        SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment(expts_to_remove{i})=[];
    end
end

NumBirds_OneTarg=length(SeqDepPitch_AcrossBirds_ONETARG.birds);


%% FILTER ALL SYLLABLES (i.e. giving each one vector of features)
% THIS CURRENTLY ONLY WORKS FOR SINGLE TARGET EXPERIMENTS

% format as a structure.

% first get all syllable names without redundancy - use motifs that I have
% hand designated


for i=1:NumBirds_OneTarg;
    NumExperiments=length(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment);
    
    for ii=1:NumExperiments;
        
        % - what motifs are in this expt?
        motifs=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder;
        
        % -- GET INFO ABOUT TARGET FOR THIS EXPT -------------------------
        targsyl=[];
        targsyl_pre=[];
        targsyl_post=[];
        % targ syl
        targsyl=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        % - preceding syl
        try % try: sometimes targyl is just one letter. sometimes is >1 letter, but is 1st in the string. that gives error
            if length(targsyl)>1;
                tmp=regexp(targsyl,'[A-Z]'); % find caps
                targsyl_pre=targsyl(tmp-1);
            end
        catch err
            disp('error: targsyl is not multiletter string, or is at start of string'); % if this is triggered, use method used for post syl, below.
            keyboard
        end
        % - post syl
        for k=1:length(motifs);
            if any(strcmp(motifs{k},targsyl));
                % find index of target
                ind=find(strcmp(motifs{k},targsyl));
                
                % post syl is one ahead
                if ind<length(motifs{k})
                    targsyl_post=motifs{k}{ind+1};
                    if length(targsyl_post)>1; % just get the single syl
                        ind=regexp(targsyl_post,'[A-Z]');
                        targsyl_post=lower(targsyl_post(ind));
                    else
                        targsyl_post=lower(targsyl_post);
                    end
                else
                    % then this is the last syl in motif, no post syl
                    targsyl_post=nan;
                end
            end
        end
        % - convert targ syl to single syllable
        if length(targsyl)>1
            % where is the upper case?
            upperind=regexp(targsyl,'[A-Z]');
            single_syl_targ=lower(targsyl(upperind));
        else
            % then could be either upper or lower
            single_syl_targ=lower(targsyl);
        end
        % - check that we got everything.
        if isempty(targsyl) | isempty(targsyl_pre) | isempty(targsyl_post);
            disp('problem, see code');
            keyboard
        end
        % --------------------------------------
        
        % -- EXTRACT INFO for each syl.
        for iii=1:length(motifs);
            syls_in_motif=motifs{iii};
            
            for j=1:length(syls_in_motif);
                
                % -- extract individual syllables
                syl=syls_in_motif{j};
                
                % what is the syl (ignoring context)?
                if length(syl)>1
                    % where is the upper case?
                    upperind=regexp(syl,'[A-Z]');
                    single_syl=lower(syl(upperind));
                else
                    % then could be either upper or lower
                    single_syl=lower(syl);
                end
                
                
                % is this the target?
                if strcmp(targsyl,syl)==1;
                    is_target=1;
                else
                    is_target=0;
                end
                
                
                % is that syl similar or different to target?
                if strcmp(single_syl,single_syl_targ)==1;
                    similar_to_targ=1;
                else
                    similar_to_targ=0;
                end
                
                
                % ID the motif
                motif_num=iii;
                
                % what is preceding syl?
                if j==1;
                    % then this is the first syl
                    preceding_syl=nan;
                else
                    preceding_syl=syls_in_motif{j-1};
                end
                
                
                % is the preceding syl similar to preceding syl of target?
                if length(preceding_syl)>1;
                    tmp=regexp(preceding_syl,'[A-Z]');
                    pre_syl=preceding_syl(tmp);
                else
                    pre_syl=preceding_syl;
                end
                
                if strcmpi(pre_syl,targsyl_pre); % if they are the same
                    presyl_similar_to_targ_presyl=1; % 1 means is similar
                else
                    presyl_similar_to_targ_presyl=0;
                end
                
                
                
                % what is the post syl?
                if j==length(syls_in_motif);
                    % then this is last syl, no post syl
                    post_syl=nan;
                else
                    post_syl=syls_in_motif{j+1};
                end
                
                % is post syl similar to target post syl?
                if length(post_syl)>1;
                    tmp=regexp(post_syl,'[A-Z]');
                    post_syl=lower(post_syl(tmp));
                else
                    post_syl=lower(post_syl);
                end
                if strcmpi(post_syl,targsyl_post); % if they are the same
                    postsyl_similar_to_targ_postsyl=1; % 1 means is similar
                else
                    postsyl_similar_to_targ_postsyl=0;
                end

                
                % if this is in same motif as targsyl, how many renditions
                % away is it? (i.e. -1 is preceding, +1 is post)
                if any(strcmp(syls_in_motif,targsyl));
                    % then this motif contains targ syl
                    targsyl_ind=find(strcmp(syls_in_motif,targsyl));
                    distance_from_targ=j-targsyl_ind;
                else
                    distance_from_targ=nan;
                end
                
                
                % -- Put all those features into the original structure
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl=single_syl;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target=is_target;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ=similar_to_targ;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).motif_num=motif_num;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl=preceding_syl;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).post_syl=post_syl;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ=distance_from_targ;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl=presyl_similar_to_targ_presyl;
                SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).postsyl_similar_to_targ_postsyl=postsyl_similar_to_targ_postsyl;
                
            end
        end
    end
end




%% GET INFORMATION ABOUT TARGET (i.e. learning at target)


for i=1:NumBirds_OneTarg;
    NumExperiments=length(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment);
    
    for ii=1:NumExperiments;
        targsyl=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};

        % WN start (1st 3 days)
             ind=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).FirstWNInd;
            Y=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).meanFF_minusBaseline(ind);
       
        % SLOT back into structure
        SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Learning_by_target.WN_start_daybins=Y;

        % WN end
        if SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.NumTargs==1;
            % only one target throughout
            ind=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).LastWNInd;
            Y=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).meanFF_minusBaseline(ind);
        elseif SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.NumTargs==0;
            % one targ epoch, then 2;
            if isfield(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
                tmpfield=fieldnames(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
                Y=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(targsyl).meanFF_minusBaseline;
            end
        end
        
        % SLOT back into structure
        SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Learning_by_target.WN_end_daybins=Y;
        
        
    end
end


%% PLOT LEARNING, SUMMARY PLOTS FOR EACH EXPERIMENT
DayBinSize=2; % days at start and end of consolid
[SeqDepPitch_AcrossBirds_ONETARG]=lt_seq_dep_pitch_ACROSSBIRDS_SummaryPlot(SeqDepPitch_AcrossBirds_ONETARG, DayBinSize);

[SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_SummaryPlot(SeqDepPitch_AcrossBirds, DayBinSize);


%% PLOT LEARNING SORTING SYLLABLES (filter based on one category only)
% this works only for categorical things (including multicategory, like
% syllable distance)

% ==== CHOOSE FILTER
Params.FilterPlot.sylID_filter{1}='similar_to_targ'; % i.e. the name of the category to filter by
Params.FilterPlot.sylID_filter{2}={0, 1}; % the possible classes of the category

% filter based on preceding syllable similarity to target
% Params.FilterPlot.sylID_filter{1}='presyl_similar_to_targ_presyl';
% Params.FilterPlot.sylID_filter{2}={0,1}; %

% filter based on post syllable similarity to target post
% Params.FilterPlot.sylID_filter{1}='postsyl_similar_to_targ_postsyl';
% Params.FilterPlot.sylID_filter{2}={0,1}; %

% ==== ADDITIONAL FILTER
Params.FilterPlot.extra_filter=0; % if 1, then filter similar vs diff
Params.FilterPlot.similar_only=0; % if 1, takes only similar, if 0, takes only different; only applies if extra_filter is on. 

% ===== What days to plot for across days analysis?
Params.FilterPlot.FirstDayToPlot=1 ; % lock start to what WN day? (1st WN day would be 1);
Params.FilterPlot.LastDayToPlot=[]; % lock end to what WN day? (5 would be 5th WN day, 'end' would be end of WN, [] would be openended);


[FILTERED_DATA, SeqDepPitch_AcrossBirds_ONETARG]=lt_seq_dep_pitch_ACROSSBIRDS_FilterPlot(Params, SeqDepPitch_AcrossBirds_ONETARG);

    


%% =========================== STORING OLD STUFF - BELOW HERE:

% %% === OBSOLETE, MOVED INTO FUNCTION ABOVE - 
% % PLOT LEARNING SORTING SYLLABLES (filter based on one category only)
% 
% % Plot what data?
% data_to_plot='WN_start_daybins';
% data_to_plot='WN_end_daybins';
% 
% 
% % sylID_filter{1}='similar_to_targ'; % i.e. the name of the category to filter by
% % sylID_filter{2}={0, 1}; % the possible classes of the category
% % this works only for categorical things (including multicategory, like
% % syllable distance)
% 
% % filter based on preceding syllable similarity to target
% sylID_filter{1}='preceding_syl';
% sylID_filter{2}={'diff','similar'}; %
% 
% % other things.
% skip_target=1; % if 1, skips collecting data for target. if 0 then collects
% skip_exceptions=1; % if 1, skips exceptions entered by hand.
% 
% FILTERED_DATA=[];
% FILTERED_DATA.filter=sylID_filter{1};
% 
% 
% for i=1:NumBirds_OneTarg;
%     NumExperiments=length(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment);
%     for ii=1:NumExperiments;
%         syl_list=fieldnames(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions);
%         
%         % -- GET INFORMATION about the target syl
%         targsyl=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
%         if length(targsyl)>1;
%             tmp=regexp(targsyl,'[A-Z]'); % find caps
%             targsyl_pre=targsyl(tmp-1);
%         end
%         
%         
%         for iii=1:length(syl_list);
%             syl=syl_list{iii};
%             
%             % -- PREPROCESSING ----------------------------------------
%             % --- Skip if this is the target syl
%             if skip_target==1;
%                 if SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
%                     continue
%                 end
%             end
%             
%             % -- Skip certain ad hoc exceptions ------------------------
%             if skip_exceptions==1;
%                 if strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.birdname,'pu11wh87') && strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.ExptID,'SyntaxDepPitchShift_cbDOWN');
%                     if strcmp(syl,'dccB'); % skip because is also targeted by WN, but not designated as targsyl.
%                         continue
%                     end
%                 end
%                 
%                 if strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.birdname,'pu11wh87') && strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
%                     if strcmp(syl,'bccB'); % skip because is also targeted by WN, but not designated as targsyl.
%                         continue
%                     end
%                 end
%             end
%             
%             if strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.birdname,'pu37wh20') && strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.ExptID,'SeqDepPitchShift');
%                 if strcmp(syl,'bccB'); % skip because is also targeted by WN, but not designated as targsyl.
%                     continue
%                 end
%             end
%             
%             
%             if strcmp(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.birdname,'pu64bk13');
%                 if strcmp(syl,'jB') || strcmp(syl,'jbB') || strcmp(syl,'jbbB') ; % skip because some get 100%Wn, some non-contingent WN (could not avoid).
%                     continue
%                 end
%             end
%       
%             
%         
%         
%         
%             % -----------------------------------------------
%             
%             % --- WHAT class is this syl?
%             sylID=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).(sylID_filter{1});
%             
%             % --- WHAT group is this data? ------------------------------
%             % determining group is slightly different for different
%             % dimensions
%             group_num=nan;
%             
%             switch sylID_filter{1};
%                 case 'preceding_syl'; % then asking whether preceding syl is similar to preceding syl of target syl
%                     % what is preceding syl?
%                     if length(sylID)>1;
%                         tmp=regexp(sylID,'[A-Z]');
%                         pre_syl=sylID(tmp);
%                         
%                     else
%                         pre_syl=sylID;
%                     end
%                     
%                     if strcmpi(pre_syl,targsyl_pre); % if they are the same
%                         group_num=2; % 2 means is similar
%                     else
%                         group_num=1;
%                     end
%                     
%                 otherwise
%                     for j=1:length(sylID_filter{2}); % use for loop because not sure if is string or numerical
%                         if sylID_filter{2}{j}==sylID;
%                             
%                             % then is in group j
%                             group_num=j;
%                             break
%                         end
%                     end
%             end
%             
%             % warn user if did not find group
%             if isnan(group_num);
%                 disp(['Warning -did not find group for bird ' num2str(i) '; experiment: ' num2str(ii) '; SKIPPING']);
%                 continue
%             end
%             
%             
%             % --- COLLECT data and group with others in this groupnum ---
%             % baseline meanFF
%             Baseline_FF_mean=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
%             
%             % 1) 3-day average to start and end WN
%             WN_Start_Ind_Binned=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).FirstWNInd; % start of WN
%             WN_start_daybins=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).rawFF{WN_Start_Ind_Binned}; % ff vals, last WN bin
%             
%             % 2) 3-day average to end WN
%             if SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.NumTargs==1; % if there was only one targ, then last WN ind is straightforward
%                 WN_End_Ind_Binned=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
%                 WNend_daybins=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).rawFF{WN_End_Ind_Binned}; % ff vals,  1st WN bin
%                 
%             elseif SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.NumTargs==0; % if started with one targ, then added another, then want last days of 1-targ, not the end of WN.
%                 if isfield(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot'); % inds should be here
%                     snapshotfield=fieldnames(SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                     
%                     WNend_daybins=SeqDepPitch_AcrossBirds_ONETARG.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(snapshotfield{1}).(syl).FFvals;
%                 else
%                     continue
%                     disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have - skipping'])
%                 end
%             end
%             
%             % 3) Get day means for all days
%             
%             % 4) ...
%             
%             
%             % -- PUT ALL DATA INTO OUTPUT STRUCTURE ---------------
%             FILTERED_DATA.groups{group_num}.groupID=sylID_filter{2}{group_num};
%             
%             if isfield(FILTERED_DATA.groups{group_num},'syllable');
%                 ind=length(FILTERED_DATA.groups{group_num}.syllable);
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.Baseline_FF_mean=Baseline_FF_mean;
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_start_daybins=WN_start_daybins;
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.WN_end_daybins=WNend_daybins;
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.syl_name=syl;
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.bird_num=i;
%                 FILTERED_DATA.groups{group_num}.syllable{ind+1}.expt_num=ii;
%             else
%                 FILTERED_DATA.groups{group_num}.syllable{1}.Baseline_FF_mean=Baseline_FF_mean;
%                 FILTERED_DATA.groups{group_num}.syllable{1}.WN_start_daybins=WN_start_daybins;
%                 FILTERED_DATA.groups{group_num}.syllable{1}.WN_end_daybins=WNend_daybins;
%                 FILTERED_DATA.groups{group_num}.syllable{1}.syl_name=syl;
%                 FILTERED_DATA.groups{group_num}.syllable{1}.bird_num=i;
%                 FILTERED_DATA.groups{group_num}.syllable{1}.expt_num=ii;
%             end
%         end
%     end
% end
% 
% 
% 
% % --- PLOT --
% % 1) Plot final learning (each syl (per experiment) contributes one datapoint);
% lt_figure; hold on;
% title(['Syls sorted by: ' sylID_filter{1} '; datapoint = one syl per expt']);
% 
% NumGroups=length(FILTERED_DATA.groups);
% Y_tot=cell(NumGroups,1); % collect data to plot
% SylNames_tot=cell(NumGroups,1);
% TargNames_tot=cell(NumGroups,1);
% 
% 
% 
% for i=1:NumGroups;
%     NumSylDataPts=length(FILTERED_DATA.groups{i}.syllable);
%     for ii=1:NumSylDataPts;
%         baseline_FF=FILTERED_DATA.groups{i}.syllable{ii}.Baseline_FF_mean;
%         FFvals=FILTERED_DATA.groups{i}.syllable{ii}.WN_end_daybins;
%         
%         FF_mean_minusbase=mean(FFvals)-baseline_FF;
%         
%         % compare to target
%         birdnum=FILTERED_DATA.groups{i}.syllable{ii}.bird_num;
%         exptnum=FILTERED_DATA.groups{i}.syllable{ii}.expt_num;
%         sylname=FILTERED_DATA.groups{i}.syllable{ii}.syl_name;
% 
%         FF_targ=SeqDepPitch_AcrossBirds_ONETARG.birds{birdnum}.experiment{exptnum}.Learning_by_target.WN_end_daybins;
%         
%         Learning_rel_targ=FF_mean_minusbase/FF_targ; % learning as fraction of learning by target (minus respective baselines)
%         
%         %         % troubleshooting
%         %         if Learning_rel_targ>0.8;
%         %             keyboard
%         %         end
%         
%         % -- COLLECT DATA
%         if SeqDepPitch_AcrossBirds_ONETARG.birds{birdnum}.experiment{exptnum}.Syl_ID_Dimensions.(sylname).similar_to_targ==1;
%             Y_tot{i}=[Y_tot{i} Learning_rel_targ];
%             
%             % collect syl name for labeling purposes
%             ind=length(SylNames_tot{i});
%             SylNames_tot{i}{ind+1}=sylname;
%             
%             % collect target name
%             targname=SeqDepPitch_AcrossBirds_ONETARG.birds{birdnum}.experiment{exptnum}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
%             TargNames_tot{i}{ind+1}=targname;
%         end
%         
%     end
%     
%     % Plot all data
%     lt_plot(i-0.2+0.4*rand(length(Y_tot{i}),1),Y_tot{i});
%     
%     % Plot syl names
%     if (1);
%         for j=1:length(Y_tot{i});
%             text(i,Y_tot{i}(j),[SylNames_tot{i}{j} '  ' TargNames_tot{i}{j}]);
%         end
%     end
%     
%     % plot mean for this group
%     Y_tot_mean=mean(Y_tot{i});
%     Y_tot_sem=lt_sem(Y_tot{i});
%     
%     errorbar(i,Y_tot_mean,Y_tot_sem,'s','MarkerFaceColor','b','MarkerSize',12);
%     
%     lt_plot_zeroline;
%     
% end
% 
% set(gca,'XTick',1:NumGroups);
% 
% 
% % -- PUT X LABELS
% switch sylID_filter{1}
%     case 'similar_to_targ'
%         set(gca,'XTickLabel',{'Different','Similar'});
%     case 'preceding_syl'
%         set(gca,'XTickLabel',sylID_filter{2});
% end
% 
% % Perform t-test
% if length(Y_tot)==2; % perform if there are two groups
%     [h,p]=ttest2(Y_tot{1},Y_tot{2});
%     % plot p-value
%     Ylim=ylim;
%     text(1.5,Ylim(2)-0.1,['p=' num2str(p)],'FontSize',12,'FontWeight','bold');
%     if p<0.05;
%         plot(1.5,Ylim(2)-0.15,'*','MarkerSize',7,'Color','r');
%     end
% end
%     
%     
%     
% 
% 
% 
% %% OLD METHOD ---- PLOT LEARNING - One plot for each experiment
% disp(' ');
% disp('NOTE: ');
% disp('For days with mult epochs, assuming 1st is 1 targ only, and DaysForSnapshot{1} is that 3-day bin');
% disp('Skipping expts with multiple targets');
% disp(' ');
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%         exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         % skip if mult targets
%         if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs>1;
%             continue
%         end
%         
%         % For this expt, plot
%         figure; hold on;
%         title(['Bird: ' birdname ', experiment: ' exptID])
%         
%         % -- SIMILAR SYLS
%         sylslist1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsSame;
%         c=1;
%         sylstoplot={};
%         for iii=1:length(sylslist1);
%             syl=sylslist1{iii};
%             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
%                 % only one target throughout
%                 ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
%                 Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).meanFF_minusBaseline(ind);
%                 Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).semFF(ind);
%             elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
%                 % one targ epoch, then 2;
%                 
%                 if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
%                     tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                     Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).meanFF_minusBaseline;
%                     Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).sem;
%                 else
%                     continue
%                     disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have -skipping'])
%                 end
%             end
%             
%             c=c+1;
%             sylstoplot=[sylstoplot syl];
%             errorbar(c,Y,Ysem,'ob');
%         end
%         
%         
%         % DIFF SYLS
%         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists,'SylsDifferent');
%             sylslist2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsDifferent;
%             for iii=1:length(sylslist2);
%                 syl=sylslist2{iii};
%                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
%                     % only one target throughout
%                     ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
%                     Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).meanFF_minusBaseline(ind);
%                     Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).semFF(ind);
%                 elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
%                     % one targ epoch, then 2;
%                     
%                     if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
%                         tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                         Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).meanFF_minusBaseline;
%                         Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).sem;
%                     else
%                         continue
%                         disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have - skipping'])
%                     end
%                 end
%                 
%                 c=c+1;
%                 sylstoplot=[sylstoplot syl];
%                 errorbar(c,Y,Ysem,'or');
%                 
%             end
%         end
%         
%         set(gca,'XTick',1:c)
%         set(gca,'XTickLabel',sylstoplot)
%         
%     end
% end
% 
% 
% %% PLOT LEARNING - SUMMARY ACROSS ALL SYLS - final learning
% disp(' ');
% disp('NOTE: ');
% disp('For days with mult epochs, assuming 1st is 1 targ only, and DaysForSnapshot{1} is that 3-day bin');
% disp('Skipping expts with multiple targets');
% disp('Assuming 1st targ syl is the one for 1st epoch');
% disp(' ');
% 
% figure; hold on;
% title('All syllables across all birds and experiments');
% 
% Y_all=struct('sim',[],'dif',[]);
% 
% for i=1:NumBirds;
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     for ii=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%         exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%         
%         % skip if mult targets
%         if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs>1;
%             continue
%         end
%         
% 
%         % -- GET INFOR ABOUT TARGET
%         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
%         
%         if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
%             % only one target throughout
%             ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).LastWNInd;
%             Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).meanFF_minusBaseline(ind);
%             Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).semFF(ind);
%         elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
%             % one targ epoch, then 2;
%             if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
%                 tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                 Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(targsyl).meanFF_minusBaseline;
%                 Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(targsyl).sem;
%             else
%                 continue
%                 disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have -skipping'])
%             end
%         end
%         
%         Ytarg=Y;
%         Ysemtarg=Ysem;
%         
% 
%         
%         % -- SIMILAR SYLS
%         sylslist1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsSame;
%         c=1;
%         sylstoplot={};
%         for iii=1:length(sylslist1);
%             syl=sylslist1{iii};
%             if strcmp(syl,targsyl);
%                 continue
%             end
%             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
%                 % only one target throughout
%                 ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
%                 Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).meanFF_minusBaseline(ind);
%                 Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).semFF(ind);
%             
%             elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
%                 % one targ epoch, then 2;
%                 if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
%                     tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                     Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).meanFF_minusBaseline;
%                     Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).sem;
%                 else
%                     continue
%                     disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have -skipping'])
%                 end
%             end
%             
% %             c=c+1;
% %             sylstoplot=[sylstoplot syl];
% %             errorbar(c,Y,Ysem,'ob');
% 
%             plot(0.7+rand*0.6,Y/Ytarg,'ob','MarkerSize',5);
%             
%             if Y/Ytarg>0.9;
%                 keyboard
%             end
%             %-- Collect Y for mean
%              Y_all.sim=[Y_all.sim Y/Ytarg];
%             
%         end
%         
%         
%         % DIFF SYLS
%         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists,'SylsDifferent');
%             sylslist2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.SylsDifferent;
%             for iii=1:length(sylslist2);
%                 syl=sylslist2{iii};
%                 
%                 if strcmp(syl,targsyl);
%                     continue
%                 end
%                 
%                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
%                     % only one target throughout
%                     ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).LastWNInd;
%                     Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).meanFF_minusBaseline(ind);
%                     Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(syl).semFF(ind);
%                 elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
%                     % one targ epoch, then 2;
%                     
%                     if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
%                         tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
%                         Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).meanFF_minusBaseline;
%                         Ysem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(syl).sem;
%                     else
%                         continue
%                         disp([num2str(i) ', ' num2str(ii) ' lacking Snapshot, but should have - skipping'])
%                     end
%                 end
%                 
% %                 c=c+1;
% %                 sylstoplot=[sylstoplot syl];
% %                 errorbar(c,Y,Ysem,'or');
%                 plot(1.7+rand*0.6,Y/Ytarg,'or','MarkerSize',5);
%                 
%             %-- Collect Y for mean
%              Y_all.dif=[Y_all.dif Y/Ytarg];
% 
%                 
%             end
%         end
%         
% %         set(gca,'XTick',1:c)
% %         set(gca,'XTickLabel',sylstoplot)
% %         
%     end
% end
% 
% % PLOT MEANS
% Y_all.simMean=mean(Y_all.sim);
% Y_all.simSEM=lt_sem(Y_all.sim);
% 
% Y_all.difMean=mean(Y_all.dif);
% Y_all.difSEM=lt_sem(Y_all.dif);
% 
% errorbar(1,Y_all.simMean,Y_all.simSEM,'ok','MarkerFaceColor','b','MarkerSize',10);
% errorbar(2,Y_all.difMean,Y_all.difSEM,'ok','MarkerFaceColor','r','MarkerSize',10);
% 
% [h,p]=ttest2(Y_all.sim,Y_all.dif)
% 
% 
% ylabel('Learning (relative to target)','FontSize',12,'FontWeight','bold');
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'Similar', 'Different'},'FontSize',12,'FontWeight','bold');
% 
% line(xlim,[0 0],'Color','k','LineStyle','--');
% 



%% SAVING
savedir='/bluejay4/lucas/across_birds/seq_dep_pitch';
currdir=pwd;

try
    cd(savedir);
catch err
    mkdir(savedir);
    cd(savedir);
end

save('SeqDepPitch_AcrossBirds', 'SeqDepPitch_AcrossBirds')
save('SeqDepPitch_AcrossBirds_ONETARG', 'SeqDepPitch_AcrossBirds_ONETARG')

cd(currdir)

