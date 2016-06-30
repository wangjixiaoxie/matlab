function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANconsolGen(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl, epochfield_input, similar_only)

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
    norm_by_targsyl=1;
end


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


%% === list names of all types of syls
if (0)
disp(' --- ');

for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1 ...
                    & SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl==0;
                
                
                disp([birdname '-' exptname '-' syl]);
            end
        end
        
        
        
    end
end
end

%% === PLOT FOR ALL EXPERIMENTS AND SYLS, TARGET REVERSION VS. EACH SYLS STATS
epochfield=epochfield_input;

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

Learning_nontarg_all=[];
X1_MUSCdividePBS_targ=[];
Y1_GenMP_MeanAllSyls=[];
X2_MUSCdividePBS_targ=[];
Y2_GenMP_AllSyls=[];
Y3_GenAFP_Mean=[];
Y4_GenAFP_All=[];
X5=[];
Y5=[];
X6=[];
Y6=[];
X7=[];
Y7=[];
X2b=[];
Y2b=[];
            X5b = [];
            Y5b= [];
Y6b=[];
X6b=[];
Y7b=[];
X7b=[];
X2c=[];
Y8=[];
X5c=[];

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
            
            % ===== STATS AT TARGET
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            FF_PBS_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs;
            FF_MUSC_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
    
            targlearndir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            FF_PBS_targ=FF_PBS_targ*targlearndir;
            FF_MUSC_targ=FF_MUSC_targ*targlearndir;           
            FF_AFP_targ=FF_PBS_targ - FF_MUSC_targ;
            
            Generalization_MP=[];
            Generalization_AFP=[];
            Generalization_Learn=[];
            AFPbias=[];
            MPbias=[];
            Learning_nontarg=[];
            for j=1:length(SylsUnique);
                    
                syl=SylsUnique{j};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    continue;
                end
                
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;

                if similar_only==0;
                    % diff only
                    if similar==1
                        continue
                    end
                elseif similar_only==1;
                    % similar only
                    if similar==0
                        continue
                    end
                elseif similar_only==2;
                    % then both
                end

% ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                FF_PBS=FF_PBS*targlearndir;
                FF_MUSC=FF_MUSC*targlearndir;
                
                % ====== CALCULATE AFP AND MUSC RELATIVE TO
                % THOSE VALUES AT TARGET (I.E. GENERALIZATIONS)
                FF_AFP=FF_PBS-FF_MUSC;
                
                Generalization_MP=[Generalization_MP FF_MUSC/FF_MUSC_targ];
                Generalization_AFP=[Generalization_AFP FF_AFP/FF_AFP_targ];
                Generalization_Learn=[Generalization_Learn FF_PBS/FF_PBS_targ];
                
                Learning_nontarg=[Learning_nontarg FF_PBS];
                AFPbias=[AFPbias FF_AFP];
                MPbias=[MPbias FF_MUSC];
            end
            
            % ==== PLOT FOR THIS EXPERIMENT
            Learning_nontarg_all=[Learning_nontarg_all Learning_nontarg];
            
            % --- 1) 
            X1_MUSCdividePBS_targ=[X1_MUSCdividePBS_targ FF_MUSC_targ/FF_PBS_targ];
            Y1_GenMP_MeanAllSyls=[Y1_GenMP_MeanAllSyls mean(Generalization_MP)];
            
                    
            % --- 2)
            X2_MUSCdividePBS_targ=[X2_MUSCdividePBS_targ FF_MUSC_targ/FF_PBS_targ*ones(1,length(Generalization_MP))];
            Y2_GenMP_AllSyls=[Y2_GenMP_AllSyls Generalization_MP];
                     
            % ---
            X2b=[X2b FF_MUSC_targ*ones(1,length(Generalization_MP))];
            Y2b=Y2_GenMP_AllSyls;
            
            X2c=[X2c FF_AFP_targ*ones(1,length(Generalization_MP))];
            
            % --- 3)
            X3=X1_MUSCdividePBS_targ;
            Y3_GenAFP_Mean=[Y3_GenAFP_Mean mean(Generalization_AFP)];
            
             % --- 4)
            X4=X2_MUSCdividePBS_targ;
            Y4_GenAFP_All=[Y4_GenAFP_All Generalization_AFP];
           
            % --- 5) 
            X5 = [X5 FF_PBS_targ* ones(1,length(Generalization_MP))];
            Y5= [Y5 AFPbias];
            
            % --- 5) 
            X5b = [X5b FF_PBS_targ*ones(1,length(Generalization_MP))];
            Y5b= [Y5b MPbias];
            
            X5c=[X5c FF_PBS_targ];
            
            % --- 6)
            X6 = [X6 FF_AFP_targ*ones(1,length(Generalization_MP))];
            Y6= Y5;
            
            % --- 6)
            X6b = [X6b FF_AFP_targ];
            Y6b= [Y6b mean(AFPbias)];

            % ---- 7) 
            X7 = [X7 FF_MUSC_targ*ones(1,length(Generalization_MP))];
            Y7= [Y7 MPbias];
            
            % ---- 7) 
            X7b = [X7b FF_MUSC_targ];
            Y7b= [Y7b mean(MPbias)];
            
            
            
            % --- 8)
            Y8=[Y8 Generalization_Learn];
            
            
            
           
%             plot(FF_MUSC_targ/FF_PBS_targ, Generalization_MP, 'bo');
%             plot(FF_MUSC_targ/FF_PBS_targ, mean(Generalization_Learn), 'o');
%             plot(FF_MUSC_targ/FF_PBS_targ, mean(Generalization_MP), 'o');
            
            
        end
    end
end


% ====== PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y1_GenMP_MeanAllSyls, X1_MUSCdividePBS_targ, 1)
ylabel('Generalization (MP)');
xlabel('Target (MUSC/PBS)');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y2_GenMP_AllSyls, X2_MUSCdividePBS_targ, 1)
ylabel('Generalization (MP)');
xlabel('Target (MUSC/PBS)');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y2b, X2b, 1)
ylabel('Generalization (MP)');
xlabel('MUSC Targ (hz)');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y3_GenAFP_Mean, X3, 1)
ylabel('Generalization (AFP)');
xlabel('Target (MUSC/PBS)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y4_GenAFP_All, X4, 1)
ylabel('Generalization (AFP)');
xlabel('Target (MUSC/PBS)');



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y8, X2c, 1)
ylabel('Generalization (Learn)');
xlabel('AFP Targ (hz)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y8, X2b, 1)
ylabel('Generalization (Learn)');
xlabel('MUSC Targ (hz)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y8, X4, 1)
ylabel('Generalization (Learn)');
xlabel('Targ (MUSC/PBS)');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Learning_nontarg_all, X2b, 1)
ylabel('Nontarg (Learn)');
xlabel('MUSC Targ (hz)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Learning_nontarg_all, X4, 1)
ylabel('Nontarg (Learn)');
xlabel('Targ (MUSC/PBS)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Learning_nontarg_all, X5, 1)
ylabel('Nontarg (Learn)');
xlabel('Targ (Learn)');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Learning_nontarg_all, X2c, 1)
ylabel('Nontarg (Learn)');
xlabel('Targ (AFP bias)');

%% ====== multilinear analysis - what is best predictor of learning at target?
% IN PROGRESS: WANT TO ESTIMATE VARIANCE EXPLAINED BY EACH PREDICTOR. THEY
% ARE CORRELATE THOUGH, SO HOW TO ESTIMATE?


XX1=X5; % targ learning
% XX2=X2c; % targ AFP bias
XX3=X2b; % targ MP bias
YY1=Learning_nontarg_all; % nontarg learning

% XXX=[ones(length(XX1),1) XX1' XX2' XX3'];
XXX=[ones(length(XX1),1) XX1' XX3'];
[b,bint,r,rint,stats]=regress(YY1', XXX);

% lt_figure; hold on;
% x1=XX1; x2=XX3; y=YY1;
% scatter3(x1,x2,y,'filled')
% hold on
% x1fit = min(x1):100:max(x1);
% x2fit = min(x2):10:max(x2);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
% mesh(X1FIT,X2FIT,YFIT)
% xlabel('Weight')
% ylabel('Horsepower')
% zlabel('MPG')
% view(50,10)


% === predictors are correlated 
lt_figure; hold on;
lt_regress(XX1, XX3, 1, 0, 1, 1, 'k');
plot(XX3, XX1, 'ok');


%% === important figs, showing that MP bias at nontarget is more related to target MP bias, than AFP
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y6, X6, 1)

ylabel('AFP bias')
xlabel('Targ AFP bias')
hsplots=[hsplots hsplot];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y6b, X6b, 1)
ylabel('AFP bias (mean over syls)')
xlabel('Targ AFP bias')
hsplots=[hsplots hsplot];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y7, X7, 1)
ylabel('MP bias')
xlabel('Targ MP bias')
hsplots=[hsplots hsplot];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y7b, X7b, 1)
ylabel('MP bias (mean over syls)')
xlabel('Targ MP bias')
hsplots=[hsplots hsplot];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y5, X5, 1)
ylabel('AFP bias')
xlabel('Raw learning (targ)')
hsplots=[hsplots hsplot];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_regress(Y5b, X5b, 1)
ylabel('MP bias')
xlabel('Raw learning (targ)')
hsplots=[hsplots hsplot];



linkaxes(hsplots, 'xy');


%% ==== compare slopes for AFP vs. AFP ; MP vs. MP
X6=X6(~isnan(X6));
Y6=Y6(~isnan(Y6));
X7=X7(~isnan(X7));
Y7=Y7(~isnan(Y7));


% === method: shuffle, each cycle perform two regressions 1) vs. Targ AFP;
% 2) vs targ MP. each DV point shuffles its actual AFP and MP into new AFP
% and MP. Each time get a difference in slope, difference in R2, etc. Use
% that for permutation test.

Yafp=Y6; % nontarg % these two vectors shuffle each time.
Ymp=Y7; % nontarg

Xafp=X6; % targ;  % these don't change
Xmp=X7;

% === get actual stats
[~,~,~,~,~,SummaryStats]=lt_regress(Yafp, Xafp, 0,0); % AFP
slopeAFP=SummaryStats.slope;
R2afp=SummaryStats.R2;
% MP
[~,~,~,~,~,SummaryStats]=lt_regress(Ymp, Xmp, 0,0); % AFP
slopeMP=SummaryStats.slope;
R2mp=SummaryStats.R2;

% ==== differences
Slope_MPminusAFP=slopeMP-slopeAFP;
R2_MPminusAFP=R2mp-R2afp;

% ===
Ncycles=5000;
inds=1:length(Yafp);
Ycombined=[Yafp' Ymp'];

Slope_MPminusAFP_perm=[];
R2_MPminusAFP_perm=[];

for i=1:Ncycles;
    
    % ==== randomly flip inds
    inds=randi(2, length(Yafp),1);
    
    Yafp_perm=[];
    Ymp_perm=[];
    for j=1:length(Xafp);
        Yafp_perm(j)=Ycombined(j, inds(j));
        Ymp_perm(j)=Ycombined(j, 3-inds(j));
    end

    % ==== perm stats
    [~,~,~,~,~,SummaryStats]=lt_regress(Yafp_perm, Xafp, 0,0); % AFP
    slopeAFP=SummaryStats.slope;
    R2afp=SummaryStats.R2;
    % MP
    [~,~,~,~,~,SummaryStats]=lt_regress(Ymp_perm, Xmp, 0,0); % AFP
    slopeMP=SummaryStats.slope;
    R2mp=SummaryStats.R2;

    % ==== differences
    Slope_MPminusAFP_perm=[Slope_MPminusAFP_perm slopeMP-slopeAFP];
    R2_MPminusAFP_perm=[R2_MPminusAFP_perm R2mp-R2afp];
end

lt_figure; hold on;
lt_subplot(3,1,1); hold on;
title('slope (MP minus AFP)');
lt_plot_histogram(Slope_MPminusAFP_perm)
line([Slope_MPminusAFP Slope_MPminusAFP], ylim, 'Color','r');

p=sum(Slope_MPminusAFP_perm>Slope_MPminusAFP)/length(Slope_MPminusAFP_perm);
lt_plot_annotation(1, ['p=' num2str(p)], 'r');


lt_subplot(3,1,2); hold on;
title('r2 (MP minus AFP)');
lt_plot_histogram(R2_MPminusAFP_perm)
line([R2_MPminusAFP R2_MPminusAFP], ylim, 'Color','r');

p=sum(R2_MPminusAFP_perm>R2_MPminusAFP)/length(R2_MPminusAFP_perm);
lt_plot_annotation(1, ['p=' num2str(p)], 'r');




%% === compare slopes for AFP vs AFP; MP vs MP
% STOPPED - below is for multilinear regresssion

% statsMP=regstats(Y7, X7);
% statsAFP=regstats(Y6, X6);
% 
% BETA=[statsMP.beta(2), statsAFP.beta(2)];
% BETACOV=[
% 
%    [p,F] = linhyptest(s.beta, s.covb, 0, [0 0 1 -1], s.tstat.dfe)
%    % result is F=4.34, p=0.0398

%% ============ AFP vs. AFP; MP vs. MP; orthogonal linear regression (takes into account error in x)
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots=[];


% --- AFP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('AFP (orthogonal regression)');
X6=X6(~isnan(X6));
Y6=Y6(~isnan(Y6));
plot(X6, Y6, 'ok');
p=linortfit2(X6, Y6);
Xlim=xlim; 
line(Xlim, [Xlim(1)*p(1)+p(2) Xlim(2)*p(1)+p(2)],'Color', 'b');
ylabel('AFP bias')
xlabel('Targ AFP bias')
lt_plot_annotation(1, ['slope=' num2str(p(1)) '; int=' num2str(p(2))], 'b');
hsplots=[hsplots hsplot];


% ---- MP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('MP (orthogonal regression)');
X7=X7(~isnan(X7));
Y7=Y7(~isnan(Y7));

plot(X7, Y7, 'ok');
p=linortfit2(X7, Y7);
Xlim=xlim; 
line(Xlim, [Xlim(1)*p(1)+p(2) Xlim(2)*p(1)+p(2)],'Color', 'b');
ylabel('MP bias')
xlabel('Targ MP bias')
lt_plot_annotation(1, ['slope=' num2str(p(1)) '; int=' num2str(p(2))], 'b');
hsplots=[hsplots hsplot];

linkaxes(hsplots, 'xy')



%% [ COLLECT DATA] - PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
% lt_figure; 
% hold on;
% title('
% line([0 0.05], [0 0.05]);

epochfield=epochfield_input;

LearningPBS_all=[];

MPbias_all=[];
AFPbias_all=[];
SimDiff_all=[];
TargStatus_all=[];
PreSylSimilar_all=[];
Expt_count_all=[];
Yexpt_all={};
Ybird_all={};
Y_PosRelTarg_All=[];

Generalization_MP_all=[];
Generalization_AFP_all=[];
Generalization_Learn_all=[];

cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];

CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];

cvPBS_alldays_ALLEXPTS={};
cvMUSC_alldays_ALLEXPTS={};


expt_count=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
        disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
            
            Y_FFmean_pbs=[];
            Y_FFmean_musc=[];
            Y_FFsem_pbs=[];
            Y_FFsem_musc=[];
            Y_syls={};
            Y_similar_diff=[];
            Y_istarg=[];
            Y_AFP_bias=[];
            Y_AcousticDist=[];
            Y_Corr=[];
            Y_presimilar=[];
            Yexpt={};
            Ybird={};
            Y_PosRelTarg=[];
            
            Y_Generalization_MP=[];
            Y_Generalization_AFP=[];
             Y_Generalization_Learn=[];
           
            
            % -- for CV stuff
            cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[];
            cvRatio_pvalue_UsingAllVals_ALLSYLS=[];
            
            CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[];
            CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[];
            
            cvPBS_alldays_ALLSYLS={};
            cvMUSC_alldays_ALLSYLS={};
            
            % ===== STATS AT TARGET
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            FF_PBS_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs;
            FF_MUSC_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
            FF_AFP_targ=FF_PBS_targ-FF_MUSC_targ;
            
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                
                % ====== CALCULATE AFP AND MUSC RELATIVE TO
                % THOSE VALUES AT TARGET (I.E. GENERALIZATIONS)
                FF_AFP=FF_PBS-FF_MUSC;
                
                Generalization_MP=FF_MUSC/FF_MUSC_targ;
                Generalization_AFP=FF_AFP/FF_AFP_targ;
                Generalization_Learn=FF_PBS/FF_PBS_targ;
                
                
                
                % ======================== calculate CV PBS and MUSC (for each day of
                % inactivattion, and mean across days)
                TvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_PBS;
                TvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_MUSC;
                
                FFvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_PBS_vals;
                FFvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_MUSC_vals;
                % convert FFvals to actual values, not diff from base
                FFvalsPBS=FFvalsPBS+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                FFvalsMUSC=FFvalsMUSC+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                
                
                % make sure tvals correspond to ffvals
                assert(length(TvalsPBS)==length(FFvalsPBS) & length(TvalsMUSC)==length(FFvalsMUSC), 'PRoblem - tvals and ffvals dont match');
                
                
                % -- for each day, calculate CV
                TvalsPBS_round=floor(TvalsPBS);
                TvalsMUSC_round=floor(TvalsMUSC);
                
                days=unique(TvalsPBS_round);
                
                cvPBS_alldays=[];
                cvMUSC_alldays=[];
                ffvals_DivideDayMean_AllDays_PBS=[];
                ffvals_DivideDayMean_AllDays_MUSC=[];
                % ffvals_MinusDayMean_AllDays_PBS=[];
                % ffvals_MinusDayMean_AllDays_MUSC=[];
                
                if isempty(days)
                    continue
                end
                
                for k=days
                    % for each day, get PBS and MUSC cv
                    
                    % --- PBS
                    inds=TvalsPBS_round==k; % only this day
                    ffvals=FFvalsPBS(inds);
                    cvPBS=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_PBS=[ffvals_MinusDayMean_AllDays_PBS ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_PBS=[ffvals_DivideDayMean_AllDays_PBS ffvals/mean(ffvals)];
                    
                    
                    
                    % --- MUSC
                    inds=TvalsMUSC_round==k;
                    ffvals=FFvalsMUSC(inds);
                    cvMUSC=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_MUSC=[ffvals_MinusDayMean_AllDays_MUSC ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_MUSC=[ffvals_DivideDayMean_AllDays_MUSC ffvals/mean(ffvals)];
                    
                    
                    % ====== save cv vals
                    cvPBS_alldays=[cvPBS_alldays cvPBS];
                    cvMUSC_alldays=[cvMUSC_alldays cvMUSC];
                    
                end
                
                % === get 1) mean of day CVs, and 2) CV over all days
                % (detrended)
%                 MeanOfDayCVs_PBS=mean(cvPBS_alldays);
%                 MeanOfDayCVs_MUSC=mean(cvMUSC_alldays);
                
                CVofAllDays_UsingValsDividedByDayMean_PBS=std(ffvals_DivideDayMean_AllDays_PBS);
                CVofAllDays_UsingValsDividedByDayMean_MUSC=std(ffvals_DivideDayMean_AllDays_MUSC);
                
                % ratio of MUSC CV to PBS CV, and whether that is
                % significant
                cvRatio_MUSCoverPBS_usingAllVals=CVofAllDays_UsingValsDividedByDayMean_MUSC/CVofAllDays_UsingValsDividedByDayMean_PBS;
                [~, p]=vartest2(ffvals_DivideDayMean_AllDays_PBS, ffvals_DivideDayMean_AllDays_MUSC, 'tail', 'right');
                disp([birdname '-' exptname '-' syl ': ' num2str(cvRatio_MUSCoverPBS_usingAllVals), '; p=' num2str(p)]);
                
%                 plot(CVofAllDays_UsingValsDividedByDayMean_MUSC, CVofAllDays_UsingValsDividedByDayMean_PBS, 'o');
                


                % ===== OUTPUT DATA
                % --- cv related stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS cvRatio_MUSCoverPBS_usingAllVals];
                cvRatio_pvalue_UsingAllVals_ALLSYLS=[cvRatio_pvalue_UsingAllVals_ALLSYLS p];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_PBS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_MUSC];
                
                cvPBS_alldays_ALLSYLS=[cvPBS_alldays_ALLSYLS cvPBS_alldays];
                cvMUSC_alldays_ALLSYLS=[cvMUSC_alldays_ALLSYLS  cvMUSC_alldays];
                
                
                
                % --- other stuff              
                Y_PosRelTarg=[Y_PosRelTarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ];
                Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
%                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];                
%                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                
                Yexpt=[Yexpt exptname(end-3:end)];
                Ybird=[Ybird birdname(1:4)];

                Expt_count_all=[Expt_count_all expt_count];
                
                
                Y_Generalization_MP=[Y_Generalization_MP Generalization_MP];
                Y_Generalization_AFP=[Y_Generalization_AFP Generalization_AFP];
                Y_Generalization_Learn=[Y_Generalization_Learn Generalization_Learn];
                
                
            end
            
            
            % ================= Flip sign if learning at targsyl is negative
            if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
                Y_FFmean_pbs=-1.*Y_FFmean_pbs;
                Y_FFmean_musc=-1.*Y_FFmean_musc;
                Y_AFP_bias=-1.*Y_AFP_bias;
            end
            
            % ========= Normalize by targsyl if desired (PBS learning
            % by taergsyl)
            if norm_by_targsyl==1;
                learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
                
                Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                Y_AFP_bias=Y_AFP_bias./learning_by_targ;
            end
            
            
            
            % ============================ COLLECT DATA TO PLOT FOR ALL
            % EXPERIMENTS
            if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
                
                % -- cv stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS];
                cvRatio_pvalue_UsingAllVals_ALLEXPTS=[cvRatio_pvalue_UsingAllVals_ALLEXPTS cvRatio_pvalue_UsingAllVals_ALLSYLS];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS];
                
                cvPBS_alldays_ALLEXPTS=[cvPBS_alldays_ALLEXPTS cvPBS_alldays_ALLSYLS];
                cvMUSC_alldays_ALLEXPTS=[cvMUSC_alldays_ALLEXPTS  cvMUSC_alldays_ALLSYLS];

                
                % -- other stuff
                LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
                MPbias_all=[MPbias_all Y_FFmean_musc];
                AFPbias_all=[AFPbias_all Y_AFP_bias];
                SimDiff_all=[SimDiff_all Y_similar_diff];
                TargStatus_all=[TargStatus_all Y_istarg];
                PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
                
                Generalization_AFP_all=[Generalization_AFP_all Y_Generalization_AFP];
                 Generalization_MP_all=[Generalization_MP_all Y_Generalization_MP];
                 Generalization_Learn_all=[Generalization_Learn_all Y_Generalization_Learn];
               
                
                
                Yexpt_all=[Yexpt_all Yexpt];
                Ybird_all=[Ybird_all Ybird];
                Y_PosRelTarg_All=[Y_PosRelTarg_All Y_PosRelTarg];
                
                
                
                expt_count=expt_count+1;
                
            end
        end
    end
end


%% AFP VS. AFP, MP VS. MP; SLOPES DIFFERENT?

keyboard


% ========= 1) USING ALL RENDS [compare means, between AFP and MP groups);
% --- AFP
X6=X6(~isnan(X6)); % targ AFP bias
Y6=Y6(~isnan(Y6)); % nontarg AFP bias


% ---- MP
X7=X7(~isnan(X7)); % targ MP bias
Y7=Y7(~isnan(Y7)); % nontarg MP bias


Xtot=[X6, X7];
Ytot=[Y6, Y7];
AFPbiasIndicator=zeros(1, length(Xtot));
AFPbiasIndicator(1:length(X6))=1;


Xname='target (hz)';
Yname='nontarg (hz)';

Zname='AFP bias? (otherwise MPbias)';

[h, a, c, stats]=aoctool(Xtot, Ytot, AFPbiasIndicator, '', Xname, Yname, Zname, 'on', 'parallel lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.04, 'estimate', 'pmm');

disp('SUGGESTS that for the same value of bias at targ, there is less AFP bias at nontarg compared to MP bias at nontarg');
disp('That is assuming same slope for AFPvsAFP and MPvsMP. Thus real effect is likely even stronger, because the slope for AFP is almost (p~0.055) smaller than the slope for MP vs. MP')
disp('so assuming same slope actually gives more close value for them than if assume different slope');
disp('this effect holds even if use separate slopes, but then that marginal mean is taking into account the mean target AFP and MP separately');
disp('this effect does not hold if use the means for each experiment, instead of using all nontargs');

pause


% ========= 1) USING ALL RENDS
% --- AFP
X6=X6(~isnan(X6)); % targ AFP bias
Y6=Y6(~isnan(Y6)); % nontarg AFP bias


% ---- MP
X7=X7(~isnan(X7)); % targ MP bias
Y7=Y7(~isnan(Y7)); % nontarg MP bias


Xtot=[X6, X7];
Ytot=[Y6, Y7];
AFPbiasIndicator=zeros(1, length(Xtot));
AFPbiasIndicator(1:length(X6))=1;


Xname='target (hz)';
Yname='nontarg (hz)';

Zname='AFP bias? (otherwise MPbias)';

[h, a, c, stats]=aoctool(Xtot, Ytot, AFPbiasIndicator, '', Xname, Yname, Zname, 'on', 'separate lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'slope');
disp('also true that, if allow separate lines, then the slope for AFP vs AFP is almost lower than the slope for MP vs. MP');
pause

[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'pmm');
disp('also true that, if allow separate lines, then the marginal mean for AFP is lower than for MP');
pause



% ====== LEARNING AND MP BIAS
X5=X5(~isnan(X5)); % targ learning
Learning_nontarg_all=Learning_nontarg_all(~isnan(Learning_nontarg_all)); % nontarg learning


% ---- MP
X7=X7(~isnan(X7)); % targ MP bias
Y7=Y7(~isnan(Y7)); % nontarg MP bias


Xtot=[X5, X7]; % targ
Ytot=[Learning_nontarg_all, Y7]; % nontarg
LearningIndicator=zeros(1, length(Xtot));
LearningIndicator(1:length(X5))=1;


Xname='target (hz)';
Yname='nontarg (hz)';

Zname='Learning1_MPbias0';

[h, a, c, stats]=aoctool(Xtot, Ytot, LearningIndicator, '', Xname, Yname, Zname, 'on', 'separate lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'slope');


% ---- plot separate linear regressions
lt_figure; hold on;
xlabel('target shift');
ylabel('nontarg shift');
title('blue: nontarg(same); bk: targ');
% --- LEARNING;
x=X5(~isnan(X5)); % targ
y=Learning_nontarg_all(~isnan(Learning_nontarg_all)); % nontarg
lt_regress(y, x, 1, 0, 1, 1, 'k');

% -- MP
x=X7(~isnan(X7)); % targ
y=Y7(~isnan(Y7)); % nontarg
lt_regress(y, x, 1, 0, 1, 1, 'b');

% ------ PLOT LINEAR REGRESSIONS ON SAME PLOT, WITH ARROW SHOWING THE
% CHANGE DUE TO INACTIVATION
lt_figure; hold on;
xlabel('target shift');
ylabel('nontarg shift');
title('blue: nontarg(same); bk: targ');
% --- LEARNING;
x=X5(~isnan(X5)); % targ
y=Learning_nontarg_all(~isnan(Learning_nontarg_all)); % nontarg
lt_plot(x, y, {'Color', 'k'});

% -- MP
x=X7(~isnan(X7)); % targ
y=Y7(~isnan(Y7)); % nontarg
lt_plot(x, y, {'Color', 'r'});

% --- arrows (MUSC minus learning)
x=X7(~isnan(X7)) - X5(~isnan(X5)); % change for targ
y=Y7(~isnan(Y7)) - Learning_nontarg_all(~isnan(Learning_nontarg_all));
% start point for each vector
xstart=X5(~isnan(X5));
ystart=Learning_nontarg_all(~isnan(Learning_nontarg_all));

% - plot lines that go up
inds=y>=0;
quiver(xstart(inds), ystart(inds), x(inds), y(inds), 0, 'Color','b');
% - plot lines that go down
inds=y<0;
quiver(xstart(inds), ystart(inds), x(inds), y(inds), 0, 'Color','m');
    
lt_plot_zeroline; lt_plot_zeroline_vert;    

% ========== PLOT JUST REGRESSION LINES AND POINTS
lt_figure; hold on;
xlabel('target shift');
ylabel('nontarg shift');
title('blue: nontarg(same); bk: targ');
% --- LEARNING;
x=X5(~isnan(X5)); % targ
y=Learning_nontarg_all(~isnan(Learning_nontarg_all)); % nontarg
lt_regress(y, x, 1, 0, 1, 1, 'k');

% -- MP
x=X7(~isnan(X7)); % targ
y=Y7(~isnan(Y7)); % nontarg
lt_regress(y, x, 1, 0, 1, 1, 'b');

lt_plot_zeroline; lt_plot_zeroline_vert;    
pause



% ========= 1) USING 1 MEAN FOR EACH EXPERIMENT
% --- AFP
inds6=~isnan(Y6b);
inds7=~isnan(Y7b);

Xtot=[X6b(inds6), X7b(inds7)]; % targ AFP bias, targ MP bias
Ytot=[Y6b(inds6), Y7b(inds7)]; % nontarg AFP bias, nontarg MP bias
AFPbiasIndicator=zeros(1, length(Xtot));
AFPbiasIndicator(1:length(inds6))=1;


Xname='target (hz)';
Yname='nontarg (hz)';

Zname='AFP bias? (otherwise MPbias)';

[h, a, c, stats]=aoctool(Xtot, Ytot, AFPbiasIndicator, '', Xname, Yname, Zname, 'on', 'separate lines');
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'estimate', 'slope');



%% ===== compare MP bias and Learning, one val per experiment, paired

Targ_Learn=[X5c];
Targ_MP=[X7b];

Nontarg_Learn=Y6b+Y7b;
Nontarg_MP=[Y7b];

inds=~isnan(Y7b);

Targ_Learn=Targ_Learn(inds);
Targ_MP=Targ_MP(inds);

Nontarg_Learn=Nontarg_Learn(inds);
Nontarg_MP=Nontarg_MP(inds);


lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('one val each expt, targ(l bar) - nontarg (r bar) paired');
xlabel('PBS ---- MUSC');
% --- PBS
xx=[1 2];
yy=[Targ_Learn' Nontarg_Learn'];

plot(xx+0.1, yy, 'o-k');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});

% --- MUSC
xx=[4 5];
yy=[Targ_MP' Nontarg_MP'];

plot(xx, yy, 'o-r');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});




%% ===== compare MP bias and Learning, all syls, not necessarily paired

Targ_Learn=[X5c];
Targ_MP=[X7b];

Nontarg_Learn=Learning_nontarg_all;
Nontarg_MP=[Y7];

lt_subplot(3,2,2); hold on;
title('one val each expt, targ(l bar) - nontarg (r bar) no paired');
xlabel('PBS --- MUSC');

% --- PBS (targ)
xx=[1];
yy=[Targ_Learn'];

plot(xx+0.1, yy, 'o-k');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});

% --- PBS (nontarg)
xx=[2];
yy=[Nontarg_Learn'];

plot(xx+0.1, yy, 'o-k');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});


% --- MUSC (targ)
xx=[4];
yy=[Targ_MP'];

plot(xx, yy, 'o-r');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});


% --- MUSC (nontarg)
xx=[5];
yy=[Nontarg_MP'];

plot(xx, yy, 'o-r');
lt_plot_bar(xx, mean(yy), {'Errors', lt_sem(yy)});










