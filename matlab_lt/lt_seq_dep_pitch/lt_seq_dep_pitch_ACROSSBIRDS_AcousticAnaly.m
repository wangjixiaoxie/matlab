function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_AcousticAnaly(SeqDepPitch_AcrossBirds, PARAMS, OverwriteSylIDs, rd28_AdHoc, dimension, TestIfDistanceSignDiff)
%% LT 12/15/15 - looking closely at whether sequence effects are strong even when accounting for structural differences

% === TO DO:
% Convincing that the acoustic distance measure is meaningful?  -- show can
% cluster, and show examples of syllables.

% Only keeping same-types that are in same cluster as the target syllable
% - plot new distributions of distances.

% Show that sequence still has effect

% Plot correlation, showing separate effects

% compare variance explained, either keeping the outliers or not.

% Plot example experiments where I controlled for similarity, and still
% show sequence effect


%%
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

    
bwidth=0.18;



%% ==== RECALCULATE ACOUSTIC DISTANCE USING PCA
NumDimensions=8;

for i=1:NumBirds;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        % ==== FIRST GET PCA Coords for target
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        fvPCAscore_Targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_PCAscore_mean;
        fvPCAscore_cov_Targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_PCAscore_cov;
        
        
        % ==== OTHER SYLS; exctract PCA coords, and get distance to target using
        % reduced dimensions (try all dimensions)
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if strcmp(syl, targsyl)
                continue
            end
            
            fvPCAscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            fvPCAscore_cov=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_cov;
            
            % --------- Measure euclidian, dot product, to target
            % [once for each dimension]
            for k=1:NumDimensions;
                dotProd=dot(fvPCAscore(1:k), fvPCAscore_Targ(1:k));
                euclDist=sqrt(sum((fvPCAscore(1:k)-fvPCAscore_Targ(1:k)).^2));
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(k).dotProd_RelTarg=dotProd;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(k).euclDist_RelTarg=euclDist;
            end
            
        end
    end
end

%% ============================== dot product [ignore, worse than euclid - note need to run euclid after this, as uses collected stats]

if (0)
    for dimension=1:8
    [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(['dot prod, dim: ' num2str(dimension)]);
    
    % COLLECT FOR HISTOGRAM
    DistancesAll=[];
    SimilarAll=[];
    PresimilarAll=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    continue
                end
                
                %                 dotProd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).dotProd_RelTarg;
                dotDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).dotProd_RelTarg;
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                
                % ===== COLLECT
                DistancesAll=[DistancesAll dotDist];
                SimilarAll=[SimilarAll similar];
                PresimilarAll=[PresimilarAll presimilar];
                
                disp([birdname '-' exptname '-' syl]);
                
               
                
            end
        end
    end
    
    % ==== PLOT HISTOGRAM
    % --- ALL
    [~, xcenters]=lt_plot_histogram(DistancesAll, '', 1, 0,'',1,'k');
    
    
    % -- different type
    inds=SimilarAll==0;
    color='r';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % -- similar type [diff presyl]
    inds=SimilarAll==1 & PresimilarAll==0;
    color='c';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % -- similar type [same presyl]
    inds=SimilarAll==1 & PresimilarAll==1;
    color='b';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % ---- put a vertical line indicating where minimum for diff types are
    inds=SimilarAll==0;
    x=max(DistancesAll(inds));
    line([x x], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
    
end
end



    
    


%% ========================= PLOT DISTRUBUTION OF DISTANCES, COMPARING THE DIFFERENT MEASURES

count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ============================== Euclidian
for dim=1:8
    % COLLECT FOR HISTOGRAM
    DistancesAll=[];
    SimilarAll=[];
    PresimilarAll=[];
    GeneralizationAll=[];
    CorrelationsSongAll=[];
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    continue
                end
                
                %                 dotProd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).dotProd_RelTarg;
                euclDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dim).euclDist_RelTarg;
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                
                generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
                
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'CORRELATIONS');
                corrSong=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
                else
                    corrSong=nan;
                end

                % ===== COLLECT
                DistancesAll=[DistancesAll euclDist];
                SimilarAll=[SimilarAll similar];
                PresimilarAll=[PresimilarAll presimilar];
                GeneralizationAll=[GeneralizationAll generalization];
                
                CorrelationsSongAll=[CorrelationsSongAll corrSong];
                disp([birdname '-' exptname '-' syl]);
                
                
%                 if euclDist>2 & similar==1 & presimilar==1;
%                     keyboard
%                 end
            end
        end
    end
    
    
    
    % ==== 1) PLOT HISTOGRAM [ALL CLASSES]
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(['euclidian, dim: ' num2str(dim)]);
    

    % --- ALL
    [~, xcenters]=lt_plot_histogram(DistancesAll, '', 1, 0,'',1,'k');
    
    
    % -- different type
    inds=SimilarAll==0;
    color='r';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % -- similar type [diff presyl]
    inds=SimilarAll==1 & PresimilarAll==0;
    color='c';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % -- similar type [same presyl]
    inds=SimilarAll==1 & PresimilarAll==1;
    color='b';
    
    Y=DistancesAll(inds);
    
    lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % ---- put a vertical line indicating where minimum for diff types are
    inds=SimilarAll==0;
    x=min(DistancesAll(inds));
    line([x x], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
    
    % ===== 2) HISTOGRAM (ALL SYLS)
    [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(['euclidian, dim: ' num2str(dim)]);

    % --- ALL
    [~, xcenters]=lt_plot_histogram(DistancesAll, '', 1, 0,'',1,'k');

    
    
    % ====== save min distance
    if dim==dimension;
            inds=SimilarAll==0;
        x=min(DistancesAll(inds));
        x_percentile=prctile(DistancesAll(inds), [2 5]);
        
        Params.acoustic.MinDistForDiffType.PCAeuclid_dim(dim)=x;
        Params.acoustic.MinDistForDiffType.PCAeuclid_2_5_prctiles_dim{dim}=x_percentile;
        
        
        disp(['SAVING min eucl dist for diff types using PCA dimension ' num2str(dim) ': min= '  num2str(x)]);
        
    end
    
end

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% BELOW, using single PCA dimensionality
%% === FOR EACH PAIR (SYL VS. TARG) DETERMINE WHETHER THEY ARE SIGNIFICANTLY DIFFERENT (USING RESAMPLING TECHNIQUE)

if TestIfDistanceSignDiff==1;
    Ncycles=200; % to shuffle


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            
            % --- get baseline PCA vectors (limit to desired
            % dimensionality)
            fvPCAallrends_targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_PCAscores(:, 1:dimension);            
            fvPCAallrends_syl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscores(:, 1:dimension);
            
            % ===== monte carlo permutation test
            fvMean_targ=mean(fvPCAallrends_targsyl,1);
            fvMean_syl=mean(fvPCAallrends_syl,1);
            DistActual=sqrt(sum((fvMean_targ-fvMean_syl).^2));
            
            assert(DistActual==SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg, 'ERROR!!!');
            
            DistPerm_All=[];
            Z=[fvPCAallrends_targsyl; fvPCAallrends_syl];
            SamplesX=size(fvPCAallrends_targsyl,1);
%             SamplesY=size(fvPCAallrends_syl,1);
            SamplesZ=size(Z,1);
            for n=1:Ncycles
                
                inds=randperm(SamplesZ);
                
                Xinds=inds(1:SamplesX);
                Yinds=inds(SamplesX+1:end);
                
                Xperm=Z(Xinds, :);
                Yperm=Z(Yinds, :);

                
                % calcalate distance
                XpermMean=mean(Xperm,1);
                YpermMean=mean(Yperm,1);
                
                DistPerm=sqrt(sum((XpermMean-YpermMean).^2));
                
                % -- save to distrubtion
                DistPerm_All=[DistPerm_All DistPerm];
            end
            
            % get p value
            p=sum(DistActual<DistPerm_All)/length(DistPerm_All);
            
            if p>0.01
                disp([birdname '-' exptname '- (targ)' targsyl '-' syl ' - p=' num2str(p)]);
            end
            
            % ===== OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(6).NULLVALS.EuclDist_RelTarg=DistPerm_All;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(6).NULLVALS.EuclDist_RelTarg_p=p;
            
            
        end
    end
end
end
    %% =============== CALCULATE NULL DISTRIBUTION OF DISTANCES FOR SAME-TYPE SYLLABLES
    
% === method: for each target syl, compare mean of N random syls to N other
% random syls. Get one distance there. Do that multiple times to get a mean
% distance. Plot mean distances across all target syllables to get estimate
% of idential syl acoustic distances.
% ==== NOTE: DIVIDING SAMPLES FOR TARGET IN HALF - this might be an overestimate of the magnitude of null distribution, might want to match sample sizes for all comparisons
% ==== NOTE: ACTUALLY NOT DIVIDING IN HALF - performing bootstrap instead,
% sampling with replacement

Ncycles=100; % will get average distance over this many cycles

DistMeansALL=[];
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        
        % --- get baseline PCA vectors (limit to desired
        % dimensionality)
        fvPCAallrends_targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_PCAscores(:, 1:dimension);
        
        Nsamps=size(fvPCAallrends_targsyl,1);
        
        
        % --- resample (with replacement) multiple times [BOOTSTRAP]
        DistShuffALL=[];
        for j=1:Ncycles;
            Xinds=randi(Nsamps, Nsamps, 1); %
            Yinds=randi(Nsamps, Nsamps, 1); %
            
            Xshuff=fvPCAallrends_targsyl(Xinds, :);
            Yshuff=fvPCAallrends_targsyl(Yinds, :);
            
            XshuffMean=mean(Xshuff);
            YshuffMean=mean(Yshuff);
            
            DistShuff=sqrt(sum((XshuffMean-YshuffMean).^2));
            
            DistShuffALL=[DistShuffALL DistShuff];
        end
        
        
%                 % --- resample (with replacement) multiple times [SPLIT IN
%                 % TWO]
%         DistShuffALL=[];
%         for j=1:Ncycles;
%             inds=randperm(Nsamps);
%             Xinds=inds(1:floor(Nsamps/2));
%             Yinds=inds(floor(Nsamps/2)+1:end);
% 
%             Xshuff=fvPCAallrends_targsyl(Xinds, :);
%             Yshuff=fvPCAallrends_targsyl(Yinds, :);
%             
%             XshuffMean=mean(Xshuff);
%             YshuffMean=mean(Yshuff);
%             
%             DistShuff=sqrt(sum((XshuffMean-YshuffMean).^2));
%             
%             DistShuffALL=[DistShuffALL DistShuff];
%         end
        
        % === calculate mean Distance
        DistMean=mean(DistShuffALL);
        
        % ==== OUTPUT
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).ACOUSTIC.usingPCA_dimensions(6).NULLVALS.EuclDist_VsSelf_Bootstrap=DistMean;

        
        % === collect across all experimetns
        DistMeansALL=[DistMeansALL DistMean];
        
    end
end
    
            
            
% === PLOT
lt_figure; hold on;
title('null distance shuffled for targs (one val for eaech targ)');
lt_plot_histogram(DistMeansALL);

            
            
            

%% CALCULATE ACOUSTIC DISTANCE FOR PRESYLS
% ==== NOTE: HAS MAJOR ISSUES - some syls presyl not defined (e.g. if
% presyl is j). Others presyl name does not mathc presyl single syl.

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        disp(' ');
        disp([' ------ ' birdname '-' exptname]);
        
        targPresyl='';
        sylPresyl='';
        
        % ====== SKIP IF TARGET LACKS PRESYL
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).preceding_syl);
            continue;
        end
        
        % --- find presyl
        motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).motif_num;
        MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{motifnum};
        PosInMotif=find(strcmp(MotifSyls, targsyl));
        
        if PosInMotif==1;
            disp(['PROBLEM - targ has presyl, but presyl does not have structure, skipping; targsyl = ' targsyl]);
            continue
        end
        
        targPresyl=MotifSyls{PosInMotif-1};
        
        % -- throw out if presyl does not match actual presyl
        if length(targPresyl)>1;
            targPresyl_lower=regexp(targPresyl, '[A-Z]', 'match');
            targPresyl_lower=lower(targPresyl_lower{1});
        else
            targPresyl_lower=targPresyl;
        end
        
        if ~strcmp(targPresyl_lower, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).preceding_syl);
            disp(['PROBLEM - presyl not acrurate, skipping; ' targPresyl_lower '-' targPresyl '-' targsyl]);
            
            continue
        end
        
        disp([targPresyl_lower '-' targPresyl '-' targsyl]);
        
        
        % ============== GO THRU ALL SYLS AND COMPARE DISTANCE OF PRESYL
        % RELATIVE TO TARGET PRESYL
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if strcmp(targsyl, syl);
                continue
            end
            
            
            % ==== DEFAULT OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_PreSylvsTargPreSyl=nan;
            
            
            % ===================================================
            % ---- SKIP IF LACKS PRESYL
            if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl);
                continue;
            end
            
            % --- find presyl
            motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).motif_num;
            MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{motifnum};
            PosInMotif=find(strcmp(MotifSyls, syl));
            
            if PosInMotif==1;
                disp(['PROBLEM - syl has presyl, but presyl does not have structure, skipping; syl = ' syl]);
                continue
            end
            sylPresyl=MotifSyls{PosInMotif-1};
            
            % -- throw out if presyl does not match actual presyl
            if length(sylPresyl)>1;
                sylPresyl_lower=regexp(sylPresyl, '[A-Z]', 'match');
                sylPresyl_lower=lower(sylPresyl_lower{1});
            else
                sylPresyl_lower=sylPresyl;
            end
            
            if ~strcmp(sylPresyl_lower, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl);
                disp(['PROBLEM - presyl not acrurate, skipping; ' sylPresyl_lower '-' sylPresyl '-' syl ' - (actual presyl)' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl]);
                
                continue
            end
            
            
            
            % ========================== SYL AND TARG BOTH HAVE VALID
            % PRESYLS - COLLECT DIFFERENCES
            targPresyl_fvPCA=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targPresyl).fv_baseline_PCAscore_mean(1:dimension);
            sylPresyl_fvPCA=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(sylPresyl).fv_baseline_PCAscore_mean(1:dimension);
            
            Distance=sqrt(sum((targPresyl_fvPCA-sylPresyl_fvPCA).^2));
            
            % ==== COLLECT OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_PreSylvsTargPreSyl=Distance;
        end
    end
end



%% CALCULATE ACOUSTIC DISTANCE FOR 2PRESYL [SAME CODE AS ABOVE!!]
% ==== NOTE: POTENTIALLY HAS MAJOR ISSUES - some syls presyl not defined (e.g. if
% presyl is j). Others presyl name does not mathc presyl single syl.

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        disp(' ');

        targPresyl='';
        sylPresyl='';
        
        % ====== SKIP IF TARGET LACKS 2-PRESYL
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).two_syl_back);
            disp([birdname '-' exptname '-' targsyl '(targ): lacks 2-presyl'])
            continue;
        end
        
        % --- find presyl
        motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).motif_num;
        MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{motifnum};
        PosInMotif=find(strcmp(MotifSyls, targsyl));
        
        if PosInMotif<=2;
            disp(['PROBLEM - targ has 2-presyl, but 2-presyl does not have structure, skipping; targsyl = ' targsyl]);
            continue
        end
        
        targPresyl=MotifSyls{PosInMotif-2};
        
        % -- throw out if presyl does not match actual presyl
        if length(targPresyl)>1;
            targPresyl_lower=regexp(targPresyl, '[A-Z]', 'match');
            targPresyl_lower=lower(targPresyl_lower{1});
        else
            targPresyl_lower=targPresyl;
        end
        
        if ~strcmp(targPresyl_lower, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).two_syl_back);
            disp(['PROBLEM - 2 presyl not acrurate, skipping; ' targPresyl_lower '-' targPresyl '-' targsyl]);
            
            continue
        end
        
        disp([birdname '-' exptname  '; 2-presyl extracted: ' targPresyl_lower '(2pre, lower)-' targPresyl '(2pre)-' targsyl '(targ)']);
        
        
        % ============== GO THRU ALL SYLS AND COMPARE DISTANCE OF PRESYL
        % RELATIVE TO TARGET PRESYL
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if strcmp(targsyl, syl);
                continue
            end
            
            
            % ==== DEFAULT OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_2PreSylvsTarg2PreSyl=nan;
            
            
            % ===================================================
            % ---- SKIP IF LACKS PRESYL
            if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back);
                disp('WHAT');
                continue;
            else
                disp('YES!!!');
            end
            
            % --- find presyl
            motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).motif_num;
            MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{motifnum};
            PosInMotif=find(strcmp(MotifSyls, syl));
            
            if PosInMotif<=2;
                disp([birdname '-' exptname ' PROBLEM - syl(' syl ') has 2presyl (' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back '), but 2presyl does not have structure, skipping']);
                continue
            end
            sylPresyl=MotifSyls{PosInMotif-2};
            
            % -- throw out if presyl does not match actual presyl
            if length(sylPresyl)>1;
                sylPresyl_lower=regexp(sylPresyl, '[A-Z]', 'match');
                sylPresyl_lower=lower(sylPresyl_lower{1});
            else
                sylPresyl_lower=sylPresyl;
            end
            
            if ~strcmp(sylPresyl_lower, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back);
                disp(['PROBLEM - 2presyl not acrurate, skipping; ' sylPresyl_lower '-' sylPresyl '-' syl ' - (actual presyl)' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back]);
                
                continue
            end
            
            
            
            % ========================== SYL AND TARG BOTH HAVE VALID
            % PRESYLS - COLLECT DIFFERENCES
            targPresyl_fvPCA=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targPresyl).fv_baseline_PCAscore_mean(1:dimension);
            sylPresyl_fvPCA=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(sylPresyl).fv_baseline_PCAscore_mean(1:dimension);
            
            Distance=sqrt(sum((targPresyl_fvPCA-sylPresyl_fvPCA).^2));
            
            % ==== COLLECT OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_2PreSylvsTarg2PreSyl=Distance;
        end
    end
end



%% ========================= COLLECT STATS


% ============================== Euclidian
    
    % COLLECT FOR HISTOGRAM
    DistancesAll=[];
    SimilarAll=[];
    PresimilarAll=[];
    GeneralizationAll=[];
    CorrelationsSongAll=[];
    PreSylDistancesAll=[];
    
    BirdnamesAll={};
    ExptNamesAll={};
    SylsAll={};
    
    
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    continue
                end
                
                %                 dotProd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).dotProd_RelTarg;
                euclDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg;
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                
                generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
                
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl), 'CORRELATIONS');
                corrSong=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
                else
                    corrSong=nan;
                end

                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension), 'euclDist_PreSylvsTargPreSyl');
                    PresylDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_PreSylvsTargPreSyl;
                    
                else
                    PresylDist=nan;
                end
                    
                
                % ===== COLLECT
                DistancesAll=[DistancesAll euclDist];
                SimilarAll=[SimilarAll similar];
                PresimilarAll=[PresimilarAll presimilar];
                GeneralizationAll=[GeneralizationAll generalization];
                
                CorrelationsSongAll=[CorrelationsSongAll corrSong];
                PreSylDistancesAll=[PreSylDistancesAll PresylDist];
                
                    BirdnamesAll=[BirdnamesAll birdname];
    ExptNamesAll=[ExptNamesAll exptname];
    SylsAll=[SylsAll syl];

                
                disp([birdname '-' exptname '-' syl]);
                
                
%                 if euclDist>2 & similar==1 & presimilar==1;
%                     keyboard
%                 end
            end
        end
    end
    
    
    %% == PLOT PRESYL DIST VS. SYL DIST [not 3d]
lt_figure; hold on;
plot_text=1;

   
% [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['using PCA, dim: ' num2str(dimension)]);
xlabel('Syl distance');
ylabel('Presyl distance');
grid on;

% --- DIFF [PRESIM]
inds=SimilarAll==0 & PresimilarAll==1;
color='r';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

% some presyl is nan - plot those on y=max anyway
Ymax=max(PreSylDistancesAll)+0.3;
Y(isnan(Y))=Ymax;

lt_plot(X, Y, {'Color', color});
if plot_text==1
    indstmp=find(inds);
    for i=1:length(indstmp);
        ii=indstmp(i);
    PlotText=[BirdnamesAll{ii}(1:4) '-' ExptNamesAll{ii}(end-2:end) '-' SylsAll{ii}];
        text(X(i), Y(i), PlotText);
    end
        
end

% --- DIFF [PREDIFF]
inds=SimilarAll==0 & PresimilarAll==0;
color='m';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);
% some presyl is nan - plot those on y=max anyway
Ymax=max(PreSylDistancesAll)+0.3;
Y(isnan(Y))=Ymax;

lt_plot(X, Y, {'Color', color});
if plot_text==1
    indstmp=find(inds);
    for i=1:length(indstmp);
        ii=indstmp(i);
    PlotText=[BirdnamesAll{ii}(1:4) '-' ExptNamesAll{ii}(end-2:end) '-' SylsAll{ii}];
        text(X(i), Y(i), PlotText);
    end
        
end


% ---- SIM (PRESIM)
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);
% some presyl is nan - plot those on y=max anyway
Ymax=max(PreSylDistancesAll)+0.3;
Y(isnan(Y))=Ymax;

lt_plot(X, Y, {'Color', color});
if plot_text==1
    indstmp=find(inds);
    for i=1:length(indstmp);
        ii=indstmp(i);
    PlotText=[BirdnamesAll{ii}(1:4) '-' ExptNamesAll{ii}(end-2:end) '-' SylsAll{ii}];
        text(X(i), Y(i), PlotText);
    end
        
end

% ---- SIM (PREDIFF)
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);
% some presyl is nan - plot those on y=max anyway
Ymax=max(PreSylDistancesAll)+0.3;
Y(isnan(Y))=Ymax;

lt_plot(X, Y, {'Color', color});
if plot_text==1
    indstmp=find(inds);
    for i=1:length(indstmp);
        ii=indstmp(i);
    PlotText=[BirdnamesAll{ii}(1:4) '-' ExptNamesAll{ii}(end-2:end) '-' SylsAll{ii}];
        text(X(i), Y(i), PlotText);
    end
        
end

    %% ===== PLOT 3D DISTRIBUTION OF PRESYL DIST VS. SYL DIST

lt_figure; hold on;
title(['[3d] using PCA, dim: ' num2str(dimension)]);
xlabel('Syl distance');
ylabel('Presyl distance');
zlabel('Generalization');
grid on;

% --- DIFF [PRESIM]
inds=SimilarAll==0 & PresimilarAll==1;
color='r';

Z=GeneralizationAll(inds);
Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot3(X, Y, Z, 'o', 'Color', color);

% --- DIFF [PREDIFF]
inds=SimilarAll==0 & PresimilarAll==0;
color='m';

Z=GeneralizationAll(inds);
Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot3(X, Y, Z, 'o', 'Color', color);


% ---- SIM (PRESIM)
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Z=GeneralizationAll(inds);
Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot3(X, Y, Z, 'o', 'Color', color);

% ---- SIM (PREDIFF)
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Z=GeneralizationAll(inds);
Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot3(X, Y, Z, 'o', 'Color', color);


% % ==== linear regression for each subclass
% % ---- DIFF
% inds=SimilarAll==0;
% color='r';
% 
% Y=CorrelationsSongAll(inds);
% X=DistancesAll(inds);
% 
% [b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);
% 
% disp(['DIFF: '])
% SummaryStats
% 
% 
% % ---- SAME PRESIM
% inds=SimilarAll==1 & PresimilarAll==1;
% color='b';
% 
% Y=CorrelationsSongAll(inds);
% X=DistancesAll(inds);
% 
% [b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);
% 
% disp(['SAME(presim): '])
% SummaryStats
% 
% % ---- SAME PREDIFF
% inds=SimilarAll==1 & PresimilarAll==0;
% color='c';
% 
% Y=CorrelationsSongAll(inds);
% X=DistancesAll(inds);
% 
% [b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);
% 
% disp(['SAME(prediff): '])
% SummaryStats


%% ==== PLOT HISTOGRAMS [SYL DISTANCES]

count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['syl distances']);

    % --- ALL
    [~, xcenters]=lt_plot_histogram(DistancesAll, '', 1, 0,'',1,'k');
    

    [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['syl distances, with ksd, bw = ' num2str(bwidth) ')']);

    
%     % ---- put a vertical line indicating where minimum for diff types are
%     inds=SimilarAll==0;
%     x=min(DistancesAll(inds));
%     line([x x], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
%     
% 

    % -- PLOT KSD
Y=DistancesAll;
[f, xi, bw]=ksdensity(Y, 'bandwidth', bwidth);
plot(xi, f, '.-');

    
    
%     % -- different type
%     inds=SimilarAll==0;
%     color='r';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
%     
%     % -- similar type [diff presyl]
%     inds=SimilarAll==1 & PresimilarAll==0;
%     color='c';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
%     
%     % -- similar type [same presyl]
%     inds=SimilarAll==1 & PresimilarAll==1;
%     color='b';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    

  %% ==== PLOT HISTOGRAMS [PRESYL DISTANCES]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['presyl distances, using pca ' num2str(dimension)]);

  
    % --- ALL
    Y=PreSylDistancesAll(~isnan(PreSylDistancesAll));
    [~, xcenters]=lt_plot_histogram(Y, '', 1, 0,'',1,'k');
    
    % ---- put a vertical line indicating where minimum for diff types are
    inds=SimilarAll==0;
    x=min(PreSylDistancesAll(inds));
    line([x x], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')

    
    
    % ==== KSD
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['presyl distances, with ksd, bw = ' num2str(bwidth) ')']);

    
%     % ---- put a vertical line indicating where minimum for diff types are
%     inds=SimilarAll==0;
%     x=min(DistancesAll(inds));
%     line([x x], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
%     
% 

    % -- PLOT KSD
    Y=PreSylDistancesAll(~isnan(PreSylDistancesAll));
[f, xi, bw]=ksdensity(Y, 'bandwidth', bwidth);
plot(xi, f, '.-');



  %% ==== PLOT HISTOGRAMS [SYL+ PRESYL DISTANCES]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('syl/presyl combined distances');

  
    % --- ALL
    Y1=PreSylDistancesAll(~isnan(PreSylDistancesAll));
    Y2=DistancesAll;
    Y=[Y1 Y2];
    
    [~, xcenters]=lt_plot_histogram(Y, '', 1, 0,'',1,'k');
    
    
    % ==== KSD
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['syl/presyl combined distances, with ksd, bw = ' num2str(bwidth)]);

    
    % -- PLOT KSD
[f, xi, bw]=ksdensity(Y, 'bandwidth', bwidth);
plot(xi, f, '.-');

    
        

%% ==== PLOT ACOUSTIC DIST VS CORRELATIONS


[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['using PCA, dim: ' num2str(dimension)]);
xlabel('acoustic dist');
ylabel('corr(song by song)');


% --- DIFF
inds=SimilarAll==0;
color='r';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});

% ---- SIM (PRESIM)
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});

% ---- SIM (PREDIFF)
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});


% ==== linear regression for each subclass
% ---- DIFF
inds=SimilarAll==0;
color='r';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['DIFF: '])
SummaryStats


% ---- SAME PRESIM
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['SAME(presim): '])
SummaryStats

% ---- SAME PREDIFF
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Y=CorrelationsSongAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['SAME(prediff): '])
SummaryStats



%% ===== PLOT generalization relative to distance
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['using PCA, dim: ' num2str(dimension)]);
xlabel('syl distance');
ylabel('generalization');

% --- DIFF
inds=SimilarAll==0;
color='r';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});

% ---- SIM (PRESIM)
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});

% ---- SIM (PREDIFF)
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

lt_plot(X, Y, {'Color', color});


% ==== linear regression for each subclass
% ---- DIFF
inds=SimilarAll==0;
color='r';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['DIFF: '])
SummaryStats


% ---- SAME PRESIM
inds=SimilarAll==1 & PresimilarAll==1;
color='b';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['SAME(presim): '])
SummaryStats

% ---- SAME PREDIFF
inds=SimilarAll==1 & PresimilarAll==0;
color='c';

Y=GeneralizationAll(inds);
X=DistancesAll(inds);

[b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y, X, 1, 0, 1, 0, color);

disp(['SAME(prediff): '])
SummaryStats



% =================== PLOT BARS TO OF MEANS
% --- DATA CLOSER THAN CLOSEST DIFF


% ---- REST OF THE DATA (GREATER THAN CLOSEST DIFF)


%% ====== DISPLAY NAMES OF ALL THE SAME-TYPES THAT ARE CLOSER THAN CLOSEST DIFF TYPE

disp(' ==================================== ')
disp(' NAMES OF SYLS CLOSER THAN CLOSEST DIFF TYPE');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        disp(' ---');
        disp([birdname '-' exptname '-' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl]);
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            %                 dotProd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).dotProd_RelTarg;
            euclDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg;
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            
            if         euclDist<Params.acoustic.MinDistForDiffType.PCAeuclid_dim(dimension);
                disp(['syl: ' syl '; dist: ' num2str(euclDist) '; genera: ' num2str(generalization) '; presim?: ' num2str(presimilar)]);
            end
            
            %                 if euclDist>2 & similar==1 & presimilar==1;
            %                     keyboard
            %                 end
        end
    end
end


        
    
%% ===== CUTOFF using global acoustic distance
% 1) One way to determine cutoff, for each get distribution of all pairwise
% distances across the dataset, see if those distances cluster. If they do,
% then choose the first cluster as defining "same-type syllables"

% 2) Method 2: within each experiemnt, choose those syls closer than all
% diff type syls - a problem is most experiments I did not calculate FV for
% all syls (just those with quanitifiable pitch)


SylDistCutoff=input('what value to use as cutoff for Syl Distance?');
PreSylDistCutoff=input('What value to use as cutoff for PreSyl Distance?');

PARAMS.SylClassify.SylDistCutoff=SylDistCutoff;
PARAMS.SylClassify.PreSylDistCutoff=PreSylDistCutoff;


% === REDEFINE ALL SYLS (SIMILAR AND PRESIM TO TARG), USING CUTOFF VALUE
SimilarAll_Redefined=nan(1, length(SimilarAll));
SimilarAll_Redefined(DistancesAll<SylDistCutoff)=1;
SimilarAll_Redefined(DistancesAll>=SylDistCutoff)=0;

PresimilarAll_Redefined=zeros(1, length(PresimilarAll)); % default is pre diff (e.g. if nan)
PresimilarAll_Redefined(PreSylDistancesAll<PreSylDistCutoff)=1;


% ======================= OUTPUT THE CUTOFF INFO
if OverwriteSylIDs==1;
    for i=1:NumBirds;
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            disp(' ---');
            disp([birdname '-' exptname '-' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl]);
            
            
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ; % save backup of hand labeled
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl_HandLab=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl; % save backup of hand labeled
                
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    
                    continue
                end
                
                % ==== similar?
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg<SylDistCutoff;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ=1;
                else
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ=0;
                end
                
                % === presimiklar?
                if rd28_AdHoc==1;
                    %  then use hand labeled
                    if strcmp(birdname, 'rd28pu64');
                        continue
                    end
                end
                
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension), 'euclDist_PreSylvsTargPreSyl');
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_PreSylvsTargPreSyl<PreSylDistCutoff;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl=1;
                    else
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl=0;
                    end
                else
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl=0;
                end
                
                % === 2presimiklar?
                
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension), 'euclDist_2PreSylvsTarg2PreSyl');
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_2PreSylvsTarg2PreSyl<PreSylDistCutoff;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ=1;
                    else
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ=0;
                    end
                else
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ=0;
                end
                
                
                
            end
        end
    end
end



%% ====== DISPLAY NAMES OF ALL SAME TYPES AND DIFF TYPES

disp(' ');
disp(' ==================================== ')
disp(' NAMES OF RECATEGORIZED SAME TYPE AND DIFF TYPES');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        disp(' ');
        disp(['== ' birdname '-' exptname '-' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl]);
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        % --- same type
        disp('-- SAME TYPE');
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            euclDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg;
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;

            if similar==1;
                disp(['syl: ' syl '; dist: ' num2str(euclDist) '; genera: ' num2str(generalization) '; presim?: ' num2str(presimilar)]);
            end
        end
        
         % --- diff type
        disp('-- DIFF TYPE');
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue
            end
            
            euclDist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).ACOUSTIC.usingPCA_dimensions(dimension).euclDist_RelTarg;
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            generalization=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;

            if similar==0;
                disp(['syl: ' syl '; dist: ' num2str(euclDist) '; genera: ' num2str(generalization) '; presim?: ' num2str(presimilar)]);
            end
        end
       
        
        
    end
end


        
%% =================================================== PLOTS USING REDEFINED SYL CATEGORIES

count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
    

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('FOLLOWING PLOTS: redefined syl categories');



%% ==== PLOT HISTOGRAMS [SYL DISTANCES]

count=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
    
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['syl distances, using pca ' num2str(dimension)]);

    % --- ALL
    [~, xcenters]=lt_plot_histogram(DistancesAll, '', 1, 0,'',1,'k');
    
    
%     % -- different type
%     inds=SimilarAll==0;
%     color='r';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
%     
%     % -- similar type [diff presyl]
%     inds=SimilarAll==1 & PresimilarAll==0;
%     color='c';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
%     
%     % -- similar type [same presyl]
%     inds=SimilarAll==1 & PresimilarAll==1;
%     color='b';
%     
%     Y=DistancesAll(inds);
%     
%     lt_plot_histogram(Y, xcenters, 1, 0, '', 1, color);
    
    % ---- put a vertical line indicating where minimum for diff types are

    line([SylDistCutoff SylDistCutoff], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
    

  %% ==== PLOT HISTOGRAMS [PRESYL DISTANCES]

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(['presyl distances, using pca ' num2str(dimension)]);

  
    % --- ALL
    Y=PreSylDistancesAll(~isnan(PreSylDistancesAll));
    [~, xcenters]=lt_plot_histogram(Y, '', 1, 0,'',1,'k');
    
    

    % ---- put a vertical line indicating where minimum for diff types are
    line([PreSylDistCutoff PreSylDistCutoff], ylim, 'Color','r','LineWidth',2, 'LineStyle','--')
    

    

%% generalziation

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('generalization');


% === same syl, same seq
inds=SimilarAll_Redefined==1 & PresimilarAll_Redefined==1;
x=1;
color='b';

Y=GeneralizationAll(inds);

plot(x+0.1, Y, 'o','Color',color);
lt_plot_bar(x, mean(Y), {'Errors',lt_sem(Y), 'Color',color});

% === same syl, Diff seq
inds=SimilarAll_Redefined==1 & PresimilarAll_Redefined==0;
x=2;
color='c';

Y=GeneralizationAll(inds);

plot(x+0.1, Y, 'o','Color',color);
lt_plot_bar(x, mean(Y), {'Errors',lt_sem(Y), 'Color',color});

% === DIFF SYL, SAME SEQ
inds=SimilarAll_Redefined==0 & PresimilarAll_Redefined==1;
x=3;
color='r';

Y=GeneralizationAll(inds);

plot(x+0.1, Y, 'o','Color',color);
lt_plot_bar(x, mean(Y), {'Errors',lt_sem(Y), 'Color',color});


% === DIFF SYL, DIFF SEQ
inds=SimilarAll_Redefined==0 & PresimilarAll_Redefined==0;
x=4;
color='m';

Y=GeneralizationAll(inds);

plot(x+0.1, Y, 'o','Color',color);
lt_plot_bar(x, mean(Y), {'Errors',lt_sem(Y), 'Color',color});


% ------
Xlabels={'Same type, same seq', 'Same type, diff seq', 'Diff type, same seq', 'Diff syl, diff seq'};
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 90);


%% Presyl dist vs. syl distance

[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Presyl dist vs. syl dist');
ylabel('Presyl dist');
xlabel('Syl dist');


% === same syl, same seq
inds=SimilarAll_Redefined==1 & PresimilarAll_Redefined==1;
color='b';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot(X, Y, 'o','Color',color);
plot(mean(X), mean(Y), 's', 'MarkerSize',8, 'Color',color, 'MarkerFaceColor',color);


% === same syl, Diff seq
inds=SimilarAll_Redefined==1 & PresimilarAll_Redefined==0;
color='c';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot(X, Y, 'o','Color',color);
plot(mean(X), nanmean(Y), 's', 'MarkerSize',8, 'Color',color, 'MarkerFaceColor',color);



% === DIFF SYL, SAME SEQ
inds=SimilarAll_Redefined==0 & PresimilarAll_Redefined==1;
color='r';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot(X, Y, 'o','Color',color);
plot(mean(X), mean(Y), 's', 'MarkerSize',8, 'Color',color, 'MarkerFaceColor',color);



% === DIFF SYL, DIFF SEQ
inds=SimilarAll_Redefined==0 & PresimilarAll_Redefined==0;
color='m';

Y=PreSylDistancesAll(inds);
X=DistancesAll(inds);

plot(X, Y, 'o','Color',color);
plot(mean(X), nanmean(Y), 's', 'MarkerSize',8, 'Color',color, 'MarkerFaceColor',color);


% % ------
% Xlabels={'Same type, same seq', 'Same type, diff seq', 'Diff type, same seq', 'Diff syl, diff seq'};
% set(gca, 'XTick', 1:4);
% set(gca, 'XTickLabel', Xlabels);
% rotateXLabels(gca, 90);


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Cluster in acoustic space (across all birds)





%% Plot distrubution of all same-type syllables in this space

% =================== PLOT all syllables across experiments
pcmean_ALL=[];
pcSD_all=[];

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
            %             plot3(pcmean(1), pcmean(2), pcmean(3),'o','Color',plotcol, 'MarkerFaceColor',plotcol)
            plot3(pcmean(1), pcmean(2), pcmean(3),'.','Color',[0.8 0.8 0.8], 'MarkerFaceColor',[0.8 0.8 0.8])
            
            % ==== COLLECT ACROSS ALL EXPERIMENTS
            pcmean_ALL=[pcmean_ALL; [pcmean(1) pcmean(2) pcmean(3)]];
            pcSD_all=[pcSD_all; [pcSD(1) pcSD(2) pcmean(3)]];
            
        end
        
        
        % -----------------------------------------------------
        
    end
end



% ================= OVERLAY SAME-TYPE SYLS AND TARGETS FOR EACH EXPERIMENT
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
            
            % ===== check that this is same-type or target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                plotcoltmp='k';
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1;
                plotcoltmp=plotcol;
            else
                continue
            end
            
            % ====== PLOT (mean + cov)
            pcmean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1
                % target
                plot3(pcmean(1), pcmean(2), pcmean(3),'s','Color',plotcol, 'MarkerFaceColor',plotcol)
            else
                plot3(pcmean(1), pcmean(2), pcmean(3),'o','Color',plotcoltmp, 'MarkerFaceColor',plotcoltmp)
            end
            
            
        end
        
        
        % -----------------------------------------------------
        
    end
end





