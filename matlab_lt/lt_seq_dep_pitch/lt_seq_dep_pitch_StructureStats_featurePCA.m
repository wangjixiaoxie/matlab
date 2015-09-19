function [Params, PCA_Struct]=lt_seq_dep_pitch_StructureStats_featurePCA(Params, AllDays_StructStatsStruct,saveON)
%% LT 4/28/15 - changed to save in same directory as all other functions


%% LT 2/5/15 - given feature vectors calculated with lt_seq_dep_pitch_StructureStats, performs PCA, and plots syl structure. Also looks at similarity scores


% Everything here is only looking at baseline

AllSylFields=fieldnames(AllDays_StructStatsStruct.IndivSyls);

% what days to look at
if strcmp(Params.PCA.epoch,'baseline');
DaysOfInterest=Params.SeqFilter.BaselineDays;
else
    DaysOfInterest=Params.PCA.epoch;
end


%% FIRST, FOR EACH FEATURE, look at values for all syls

% GET mean baseline feature vector
for i=1:length(AllSylFields);
    syl=AllSylFields{i};

    bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DaysOfInterest); % inds of rends that are from baseline
    fvec=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,:); % these are baseline, all renditions
    
    % mean and std
    AllFeatVecMean(i,:)=mean(fvec,1);
    AllFeatVecSD(i,:)=std(fvec);
    AllFeatVecN(i)=size(fvec,1);
end

% PLOT - one plot for each feature vector, compare syls;
for i=1:length(Params.StructureStats.FeatureLegend); % for each feature
    feature=Params.StructureStats.FeatureLegend{i};
    
    figure; hold on;
    title(['Mean (SD) of ' feature ' (during baseline)']);
    
    for ii=1:length(AllSylFields);
        syl=AllSylFields{ii};
        
        if sum(strcmp(syl,Params.SeqFilter.SylLists.TargetSyls))==1; % color based on relationship to targ syl.
            
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'ok','MarkerSize',9);
        elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsSame))==1;
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'og','MarkerSize',9);
        elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsDifferent))==1;
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'or','MarkerSize',9);
        else
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'ob','MarkerSize',9);
            disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
        end
        
    end
    
    set(gca,'XTick',1:length(AllSylFields));
    
    set(gca,'XTickLabel',AllSylFields);
    
end




%% PCA of vectors

% FIRST, compile all data into one vector
% columns: all renditions of all syllables, across baseline days
% rows: feature values

AllSingleSyls=Params.DayRawDat.syllables; % use just single syls, so will not double count syls.
syllist=AllSingleSyls;
AllSylFields=fieldnames(AllDays_StructStatsStruct.IndivSyls);

PCA.AllSylsRends_vals=[];

for i=1:length(syllist);
    syl=syllist{i};
    
    % find baseline days
    bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DaysOfInterest); % inds of rends that are from baseline

    % extract baseline data
    tmp=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,:); % all rends data (baseline)
    PCA.AllSylsRends_vals=[PCA.AllSylsRends_vals; tmp];
    
end


% CONVERT TO Z-SCORE
% for each feature, get global mean (across syls and all rends), and get
% z-score rel to that
PCA.AllSylsRends_mean=mean(PCA.AllSylsRends_vals,1);
PCA.AllSylsRends_std=std(PCA.AllSylsRends_vals);

N=size(PCA.AllSylsRends_vals,1);
PCA.AllSylsRends_zscore=(PCA.AllSylsRends_vals-repmat(PCA.AllSylsRends_mean,N,1))./repmat(PCA.AllSylsRends_std,N,1);


% PERFORM PCA
[coeff, score, latent, tsquared, explained]=pca(PCA.AllSylsRends_zscore);

PCA.outputs.coeff=coeff; % orthonormal vector coefficients
PCA.outputs.score=score; % coords of original data in new coord system
PCA.outputs.explained=explained; % explained var
PCA.outputs.latent=latent; % variance explained of each PC
PCA.outputs.tsquared=tsquared; % hotelling t distance for each rend (i.e. var weighted distance from global centroid).


% PLOT 1ST VS. 2ND PC, and 1st 3 in 3d
figure; hold on;
plot3(score(:,1), score(:,2), score(:,3),'+');
title('PC3 vs. PC2 vs. PC1');

% DISPLAY PCA stats
figure; hold on;
pareto(explained);
title('% variance explained by each PC');

% WHAT DO the PCs look like?
figure; hold on;
title('PC coefficients (lines) for first 3 PCs, and all scores from data (normalized)');
biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',Params.StructureStats.FeatureLegend)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');



%% GET COORDINATES FOR SEQ-DEP SYLS


for i=1:length(AllSylFields);
    syl=AllSylFields{i};
    
    
    % find baseline days
    bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,DaysOfInterest); % inds of rends that are from baseline
        
    % first get z-score, based on mean and std that was computed over all
    % syls (above)
    Vals= AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,:); % raw vals
    N=size(Vals,1);
    Vals_z=(Vals-repmat(PCA.AllSylsRends_mean,N,1))./repmat(PCA.AllSylsRends_std,N,1);
    
    % Get PCA scores
    PCA.IndivSyls.(syl).score=Vals_z*coeff; % LT verified that coeff should not be transposed.     
    PCA.IndivSyls.(syl).Vals=Vals;
    PCA.IndivSyls.(syl).Vals_z=Vals_z;

end


% GET MEAN AND COV for each syl
for i=1:length(AllSylFields);
    syl=AllSylFields{i};
    
    PCA.IndivSyls.(syl).score_mean=mean(PCA.IndivSyls.(syl).score,1); % mean score
    PCA.IndivSyls.(syl).score_cov=cov(PCA.IndivSyls.(syl).score);   
end


% 1) PLOT scatter all vals
figcol=[0.9 0.9 0.9];
figure; hold on;
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

% first plot all
plot3(score(:,1), score(:,2), score(:,3),'ok');
title('PC3 vs. PC2 vs. PC1');

% then plot all seq dep syls
syllist=[Params.SeqFilter.SylLists.SylsSame, Params.SeqFilter.SylLists.SylsDifferent, Params.SeqFilter.SylLists.TargetSyls]; % gets all with no overlap.
PlotColors=lt_make_plot_colors(length(syllist),0);

for i=1:length(syllist);
    syl=syllist{i};
    
    score=PCA.IndivSyls.(syl).score;
    hfig2(i)=plot3(score(:,1), score(:,2), score(:,3),'o','Color',PlotColors{i},'MarkerFaceColor',PlotColors{i});
end

legend(hfig2,syllist);
set(gca,'Color',figcol)
    

% 2) PLOT MEAN AND COV
figure; hold on;
title('Mean (SD) PC scores for each syl');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

for i=1:length(syllist);
    syl=syllist{i};
    
    % plot ellipsoid of first 3 PCs
    pcmean=PCA.IndivSyls.(syl).score_mean;
    pcSD=sqrt(diag(PCA.IndivSyls.(syl).score_cov));
    
    [x,y,z]=ellipsoid(pcmean(1),pcmean(2),pcmean(3), pcSD(1), pcSD(2), pcSD(3));
    
    
    hellip(i)=surf(x,y,z,i*ones(length(z)));

    PlotColorsMat=cell2mat(PlotColors);
    PlotColorsMat=reshape(PlotColorsMat,length(PlotColorsMat)/3,3);
    colormap(PlotColorsMat)
end

legend(hellip,syllist)


%% CALCULATE SYL SIMILARITY




%% SAVE

if saveON==1;
tstamp=lt_get_timestamp(0);

cd(Params.SeqFilter.savedir);

save('Params.mat','Params');
PCA_Struct=PCA;
save('PCA_Struct.mat','PCA_Struct');

% write a text file that tells you when files were made
fid1=fopen(['DONE_featurePCA_' tstamp '.txt'],'w');
fclose(fid1);

try 
    cd FIGURES/featurePCA
catch err
       mkdir FIGURES/featurePCA
     cd FIGURES/featurePCA
end

lt_save_all_figs;

cd ../../
end


