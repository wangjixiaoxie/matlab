function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_CORR_and_PCA(SeqDepPitch_AcrossBirds, PARAMS, plotPCA_stuff, Zscore_use_rends, load_old_vector_space)

%% LT - gets correlations and acoustic distance info, 7/16/15

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% GET CORRELATIONS and REGULAR EXPRESSIONS
% First get regular expressions extracted data as that contains motif
% classified information, which allows us to perform correlation analyses
% easily.
% FOR EACH EXPT, RUN REGEXPR OR LOAD - EITHER WAY PLOT FOR VISUALIZATION

currdir=pwd;
PARAMS.global.currdir=currdir;

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        close all;
        
        % ====================== LOAD RegExpr info (which contains
        % correlations info
        
        % 1) Check whether RegExpr are already done
        cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
        
        tmp1=dir('DONE_RegExpr*');
        tmp2=dir('DONE_PlotRegExpr*');
        
        if ~isempty(tmp2) && ~isempty(tmp1);
            % then have already processed corr data. Simply load and
            % extract
            
            % -------------- EXTRACT
            clear params;
            clear data;
            
            params=load('Params');
            data=load('AllDays_RegExpr');
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=params.Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=data.AllDays_RegExpr;
            
            disp(['Extracted already-done RegExpr data for bird ' num2str(i) ', expt ' num2str(ii)]);
            
        elseif isempty(tmp2) && ~isempty(tmp1);
            % then have extracted RegExpr data, but not collected corr data
            
            % ---------------- RUN corr and plot analysis
            clear Params
            clear AllDays_RegExpr
            load('Params');
            load('AllDays_RegExpr');
            
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr, saveON, 0, 1);
            
            % ------------- SAVE TO OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
            disp(['Performed RegExprPlot analysis for bird ' num2str(i) ', expt ' num2str(ii)]);
            
        elseif isempty(tmp2) && isempty(tmp1);
            % Then have to start from scratch
            clear Params;
            clear AllDays_RegExpr;
            clear AllDays_RawDatStruct;
            
            load('Params');
            load('AllDays_RawDatStruct');
            load('AllDays_PlotLearning');
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                DoLMAN=1;
            else
                DoLMAN=0;
            end
            
            % ------ RUN REG EXPR
            Params.RegExpr.expressions=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
            saveON=1;
            [Params, AllDays_RegExpr] = lt_seq_dep_pitch_RegExpr(Params, AllDays_RawDatStruct, saveON, DoLMAN, AllDays_PlotLearning);
            
            Params.PlotRegExpr.plotWNdays=1;
            saveON=1;
            [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr, saveON, 0, 1);
            
            
            % ---------- SAVE OUTPUTS
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params=Params;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr=AllDays_RegExpr;
            
            disp(['Extracted RegExpr data and performed RegExprPlot analysis for bird ' num2str(i) ', expt ' num2str(ii)]);
            
        end
    end
end

close all;
 
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
%     if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning, 'SylFields_Unique'); % unique syls
%         syls_unique={}; % then collect unqiue syls from scratch
%         for j=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder);
%             syls_unique=[syls_unique SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{j}];
%         end
%         SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=syls_unique;
%     end
    
    % Collect stuff
%     syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{expt_ind}.INFORMATION.SylFields_Unique;
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
        tmp=dir('DONE_StructureStats_*');
        if ~isempty(tmp);
            % THEN SAVED STRUCTURE EXISTS - Load and extract data
            load AllDays_StructStatsStruct;
            load Params;
            
        else
            % THEN HAVE TO RUN CODE TO GET STATS STRUCT
            load Params;
            load AllDays_RawDatStruct;
            load AllDays_PlotLearning;
                DoLMAN=0;
            
            % NOTE: if this is LMAN expt - this code DOES only take PBS
            % (but takes all PBS, not just in time window)
            [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct, DoLMAN, AllDays_PlotLearning);
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

    
%% ====== PERFORM PCA ON ALL REPRESENTATIVE SYLLABLES (One for each bird, to construct the PCA space)
% == First collect all syllables into one large array (subsample to get same
% number of samples for each syllable)

% 1) ============= SAMPLE SIZE 
% First, display sample sizes:
disp(' ');
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
    disp(' ');
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

if plotPCA_stuff==1;
lt_figure; hold on;
title('histogram of sample sizes over all syls (one set per bird)');
hist(N_all,50);
xlabel('sample size');
ylabel('count');
end

% Third, choose the max of 50 and lowest 5 percentile
N_subsample = floor(max([50, prctile(N_all, 5)]));


% 2) ============== Subsample and collect all data
PCA.raw_data.FV_AllSylsAllBirds=[];
PCA.raw_data.FV_AllSylsAllBirds_sylmeans=[]; % each syl contributes one vector.

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
        PCA.raw_data.FV_AllSylsAllBirds_sylmeans=[PCA.raw_data.FV_AllSylsAllBirds_sylmeans; mean(featurevect_temp)];
    end
end

PCA.raw_data.N_subsample=N_subsample;


% ==== SAVE VECTOR SPACE
% - only save if not trying to load old
    pathdir='/bluejay4/lucas/across_birds/seq_dep_pitch/ACOUSTIC_VECT_SPACES';
if load_old_vector_space==0
    cd(pathdir);
    tstamp=lt_get_timestamp(0);
    
    filenm=['PCA_struct_' tstamp];
    save(filenm, 'PCA');
end

%% =============== z-score all data (relative to the global mean + SD above
% ====== DECIDE WHETHER TO USE THE CURRENT VECTOR SPACE, OR LOAD AN OLD
% VECTOR SPACE
if load_old_vector_space==1
    % confirm
    if strcmp(input('sure you want to load old acoustic vector space? (y or n)', 's'), 'y');
        disp('loading...');
        cd(pathdir)
        
        tmpdir=dir('PCA_struct*');
        
        % show the latest 10
        if length(tmpdir)>10;
            tmpdir=tmpdir(end-9:end);
        end
        
        for j=1:length(tmpdir);
            disp([num2str(j) ' -- ' tmpdir(j).name ' -- ' datestr(tmpdir(j).datenum, 'ddmmmyyyy-HHMM')]);
        end
        
        indtoload=input('WHICH FILE TO LOAD? (e.g. 1) ');
        
        clear PCA
        load(tmpdir(j).name);
        disp(['-- LOADED: ' tmpdir(j).name]);
        pause(2);
        
    end
end

% ==== RUN 
if Zscore_use_rends==0 % defualt is 1 (uses all rends). if 0, then uses mean for each syl
    PCA.raw_data.FV_AllSylsAllBirds_allrends=PCA.raw_data.FV_AllSylsAllBirds; % save backup of all rends data
    PCA.raw_data.FV_AllSylsAllBirds=PCA.raw_data.FV_AllSylsAllBirds_sylmeans; % replace all rends field with syl means
end
PCA.raw_data.FV_mean=mean(PCA.raw_data.FV_AllSylsAllBirds,1);
PCA.raw_data.FV_std=std(PCA.raw_data.FV_AllSylsAllBirds);

N=size(PCA.raw_data.FV_AllSylsAllBirds,1);

PCA.raw_data.FV_MinusMean=PCA.raw_data.FV_AllSylsAllBirds-repmat(PCA.raw_data.FV_mean,N,1);
PCA.raw_data.FV_zscore=PCA.raw_data.FV_MinusMean./repmat(PCA.raw_data.FV_std,N,1);
    



% ================== PCA 
[coeff, score, latent, tsquared, explained]=pca(PCA.raw_data.FV_zscore);

PCA.outputs.coeff=coeff; % orthonormal vector coefficients
PCA.outputs.score=score; % coords of original data in new coord system
PCA.outputs.explained=explained; % explained var
PCA.outputs.latent=latent; % variance explained of each PC
PCA.outputs.tsquared=tsquared; % hotelling t distance for each rend (i.e. var weighted distance from global centroid).

if plotPCA_stuff==1;

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
end


%% =============== FOR EACH SYL (for all experiments), get coordinates in PC space, and also z-score space.
% ==== LOAD STRUCTURESTATS (PBS)
% in PCA space
% N_dimensions=input('How many PCA dimensions to use? '); 
% NOTE ---> THIS HAS NOT BEEN IMPLEMENTED YET!!!! (currently takes all
% dimensions, so is the same as z-score space)

for i=1:NumBirds;
%     disp(['bird: ' num2str(i)]);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    num_experiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%     lt_figure; hold on; % for this bird
    
    for ii=1:num_experiments;
%         disp(['exptID: ' num2str(ii)]);
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % =============== Extract baseline FV
        save_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir;
        
        
        % ============ COLLECT FEATURE VECTORS FOR ALL BASELINE RENDITIONS
        cd(save_dir)
        % --- Does saved structure already exist?
        % IS THE FOLDER STRUCTURED IN THE OLD OR NEW FORMAT?
        tmp=dir('AllDays_RawDatStruct.mat');
        if isempty(tmp)
            % Old format -->
            disp('PROBLEM - folder structured in old format - bad...');
        else
            % New format -->
            tmp=dir('DONE_StructureStats_*');
            if ~isempty(tmp);
                % THEN SAVED STRUCTURE EXISTS - Load and extract data
                load AllDays_StructStatsStruct;
                load Params;
                
            else
                % THEN HAVE TO RUN CODE TO GET STATS STRUCT
                load Params;
                load AllDays_RawDatStruct;
                load AllDays_PlotLearning;
                DoLMAN=0;
                [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct, DoLMAN, AllDays_PlotLearning);
                close all;
            end
        end
        cd(currdir);
        % ------------------------------------------------------------
        
        % ===================== COLLECT BASELINE FVs for all unique syls
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        baseline_days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            baselineInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd, baseline_days);
            
            fv=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(baselineInds,:);
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_all_rends=fv;
        end
        
        
        % ===== PLOT EACH SYL, ONE FOR EACH BIRD/EXPT, IN 3D PCA SPACE 
%         hellip=[];
%         PlotColors=lt_make_plot_colors(length(syls_unique), 0, 0);
%         lt_subplot(ceil(num_experiments/2), 2, ii); hold on;
%         xlabel('PC1');
%         ylabel('PC2');
%         zlabel('PC3');
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            fv_array=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_all_rends;
            
            % make sure have enough sample size
            N_samps=size(fv_array,1);
            if N_samps<10;
                disp(['WARNING: sample size for bird/expt/syl: ' num2str(i) num2str(ii) syl ' is ' num2str(size(fv_array,1))]);
            end
            
            % - get z-score for fvals
            meanArray=repmat(PCA.raw_data.FV_mean, N_samps, 1);
            stdArray=repmat(PCA.raw_data.FV_std, N_samps, 1);
            
            fv_array_zscore=(fv_array-meanArray)./stdArray;
            
            % -- dot fv array with coefficients to convert to PCA coords.
            fv_array_PCAcoords=fv_array_zscore*PCA.outputs.coeff;
            
            % ==== PUT INTO OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore=fv_array_zscore;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean=mean(fv_array_zscore,1);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscores=fv_array_PCAcoords;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean=mean(fv_array_PCAcoords,1);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov=cov(fv_array_PCAcoords,1);
            
            
            % ====== PLOT (mean + cov)
%             pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
%             pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
%             
%             [x,y,z]=ellipsoid(pcmean(1),pcmean(2),pcmean(3), pcSD(1), pcSD(2), pcSD(3));
%             
%             hellip(j)=surf(x,y,z,j*ones(length(z)));
%             
%             PlotColorsMat=cell2mat(PlotColors);
%             PlotColorsMat=reshape(PlotColorsMat,length(PlotColorsMat)/3,3);
%             colormap(PlotColorsMat)
        end
        
%         title(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID);
%         legend(hellip,syls_unique);
        % -----------------------------------------------------
        
        
        
        
    end
%     lt_subtitle(SeqDepPitch_AcrossBirds.birds{i}.birdname)
    
end



%% LOAD STRUCTURESTATS (MUSC) [SAME AS ABOVE BUT MUSC]

for i=1:NumBirds;
%     disp(['bird: ' num2str(i)]);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    num_experiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
%     lt_figure; hold on; % for this bird
    
    for ii=1:num_experiments;
%         disp(['exptID: ' num2str(ii)]);
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % === check if is LMAN expt
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
            continue
        end
        
        % =============== Extract baseline FV
        save_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir;
        
        
        % ============ COLLECT FEATURE VECTORS FOR ALL BASELINE RENDITIONS
        cd(save_dir)
        % --- Does saved structure already exist?
        % IS THE FOLDER STRUCTURED IN THE OLD OR NEW FORMAT?
        tmp=dir('AllDays_RawDatStruct.mat');
        if isempty(tmp)
            % Old format -->
            disp('PROBLEM - folder structured in old format - bad...');
        else
            % New format -->
            tmp=dir('DONE_StructureStatsMUSC*');
            if ~isempty(tmp);
                % THEN SAVED STRUCTURE EXISTS - Load and extract data
                load('AllDays_StructStatsStruct_MUSC');
                AllDays_StructStatsStruct_MUSC=AllDays_StructStatsStruct;
                load Params;
            else
                % THEN HAVE TO RUN CODE TO GET STATS STRUCT
                load Params;
                load AllDays_RawDatStruct;
                load AllDays_PlotLearning;
                DoLMAN=1;
                [Params, AllDays_StructStatsStruct_MUSC]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct, DoLMAN, AllDays_PlotLearning);
                close all;
            end
        end
        cd(currdir);
        % ------------------------------------------------------------
        
        % ===================== COLLECT BASELINE FVs for all unique syls
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        baseline_days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays;
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            baselineInds=ismember(AllDays_StructStatsStruct_MUSC.IndivSyls.(syl).AllRends.Info.dayInd, baseline_days);
            
            fv=AllDays_StructStatsStruct_MUSC.IndivSyls.(syl).AllRends.FeatVect(baselineInds,:);
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_all_rends=fv;
        end
        
        
        % ===== PLOT EACH SYL, ONE FOR EACH BIRD/EXPT, IN 3D PCA SPACE 
%         hellip=[];
%         PlotColors=lt_make_plot_colors(length(syls_unique), 0, 0);
%         lt_subplot(ceil(num_experiments/2), 2, ii); hold on;
%         xlabel('PC1');
%         ylabel('PC2');
%         zlabel('PC3');
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            fv_array=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_all_rends;
            
            % make sure have enough sample size
            N_samps=size(fv_array,1);
            if N_samps<10;
                disp(['WARNING: sample size for bird/expt/syl: ' num2str(i) num2str(ii) syl ' is ' num2str(size(fv_array,1))]);
            end
            
            % - get z-score for fvals
            meanArray=repmat(PCA.raw_data.FV_mean, N_samps, 1);
            stdArray=repmat(PCA.raw_data.FV_std, N_samps, 1);
            
            fv_array_zscore=(fv_array-meanArray)./stdArray;
            
            % -- dot fv array with coefficients to convert to PCA coords.
            fv_array_PCAcoords=fv_array_zscore*PCA.outputs.coeff;
            
            % ==== PUT INTO OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_zscore=fv_array_zscore;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_zscore_mean=mean(fv_array_zscore,1);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_PCAscores=fv_array_PCAcoords;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_PCAscore_mean=mean(fv_array_PCAcoords,1);
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_PCAscore_cov=cov(fv_array_PCAcoords,1);
            
            
            % ====== PLOT (mean + cov)
%             pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
%             pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
%             
%             [x,y,z]=ellipsoid(pcmean(1),pcmean(2),pcmean(3), pcSD(1), pcSD(2), pcSD(3));
%             
%             hellip(j)=surf(x,y,z,j*ones(length(z)));
%             
%             PlotColorsMat=cell2mat(PlotColors);
%             PlotColorsMat=reshape(PlotColorsMat,length(PlotColorsMat)/3,3);
%             colormap(PlotColorsMat)
        end
        
%         title(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID);
%         legend(hellip,syls_unique);
        % -----------------------------------------------------
        
    end
%     lt_subtitle(SeqDepPitch_AcrossBirds.birds{i}.birdname)
    
end




%% =============== PLOT ALL SYLS/ALL EXPTS, in one plot [PBS]
% in PCA space
% N_dimensions=input('How many PCA dimensions to use? '); 
% NOTE ---> THIS HAS NOT BEEN IMPLEMENTED YET!!!! (currently takes all
% dimensions, so is the same as z-score space)
if plotPCA_stuff==1;

pcmean_ALL=[];
pcSD_all=[];


% 1) PLOT ALL SYLS/ALL EXPERIMENTS
lt_figure; hold on;
title('All syls, all expts (each expt diff col)');
        xlabel('PC1');
        ylabel('PC2');
        zlabel('PC3')
        
for i=1:NumBirds;
%     disp(['bird: ' num2str(i)]);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    num_experiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);

    for ii=1:num_experiments;
%         disp(['exptID: ' num2str(ii)]);
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
       plotcol=[rand rand rand];
        % ===== PLOT EACH SYL, ONE FOR EACH BIRD/EXPT, IN 3D PCA SPACE 
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
            
            % ====== PLOT (mean + cov)
            pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
            
%             lt_plot(pcmean(1), pcmean(2), {'Color',plotcol});
            plot3(pcmean(1), pcmean(2), pcmean(3),'o','Color',plotcol, 'MarkerFaceColor',plotcol)
            
        % ==== COLLECT ACROSS ALL EXPERIMENTS
        pcmean_ALL=[pcmean_ALL; [pcmean(1) pcmean(2) pcmean(3)]];
        pcSD_all=[pcSD_all; [pcSD(1) pcSD(2) pcmean(3)]];
        
        end
        
        
        % -----------------------------------------------------
        
    end
end

end
%% =============== PLOT EACH EXPERIMENT SEPARATELY, OVERLAYED WITH ALL SYLS IN ALL EXPTS
% in PCA space
% N_dimensions=input('How many PCA dimensions to use? '); 
% NOTE ---> THIS HAS NOT BEEN IMPLEMENTED YET!!!! (currently takes all
% dimensions, so is the same as z-score space)


% 1) PLOT ALL SYLS/ALL EXPERIMENTS
% lt_figure; hold on;
% title('All syls, all expts (each diff col)');
% xlabel('PC1');
% ylabel('PC2');
if plotPCA_stuff==1;



for i=1:NumBirds;
%     disp(['bird: ' num2str(i)]);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    num_experiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    

    count=1;
    SubplotsPerFig=6;
    subplotrows=3;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
        
    for ii=1:num_experiments;
%         disp(['exptID: ' num2str(ii)]);
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title(exptname);
        xlabel('PC1');
        ylabel('PC2');
        zlabel('PC3');

        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                
        PlotCols=lt_make_plot_colors(length(syls_unique), 0, 0);

        % 1) ====== PLOT ALL SYLS ACROSS ALL EXPERIMENTS FIRST.
%         lt_plot(pcmean_ALL(:,1), pcmean_ALL(:,2), {'Color', [0.8 0.8 0.8]});
        plot3(pcmean_ALL(:,1), pcmean_ALL(:,2), pcmean_ALL(:,3),'o','Color',[0.8 0.8 0.8], 'MarkerFaceColor',[0.8 0.8 0.8])

        % 2 ====== PLOT EACH SYL FOR THIS EXPT
        hplot=[];
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
            
            % ====== PLOT (mean + cov)
            pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            pcSD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
            
%             hplot(j)=lt_plot(pcmean(1), pcmean(2), {'Color',PlotCols{j}});
            hplot(j)=plot3(pcmean(1), pcmean(2), pcmean(3), 'o','Color',PlotCols{j}, 'MarkerFaceColor',PlotCols{j});
           
            
            %             [x,y,z]=ellipsoid(pcmean(1),pcmean(2),pcmean(3), pcSD(1), pcSD(2), pcSD(3));
            %
            %             hellip(j)=surf(x,y,z,j*ones(length(z)));
            %
            %             PlotColorsMat=cell2mat(PlotColors);
            %             PlotColorsMat=reshape(PlotColorsMat,length(PlotColorsMat)/3,3);
            %             colormap(PlotColorsMat)
        end
        
        legend(hplot,syls_unique);
        % -----------------------------------------------------
        
    end
    lt_subtitle(birdname)
    
end
end



