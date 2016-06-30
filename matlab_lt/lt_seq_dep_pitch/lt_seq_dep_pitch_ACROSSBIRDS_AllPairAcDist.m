function lt_seq_dep_pitch_ACROSSBIRDS_AllPairAcDist(SeqDepPitch_AcrossBirds, Params, IncludeThrownOutSyls, threshold)


%% 

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
bwidth=0.13;

%%
AllPairwiseDist=[];
figcount=1;
subplotrows=3;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];

% Determine colormap such that all pairs below threshold are shown in one
% color
clims=[threshold 2*threshold]; % i.e. dynamic range limited to this (outside all plotted as one color)


for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        if IncludeThrownOutSyls==1;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        else
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique; 
        end
        disp(' ----' );
        disp(SylsUnique);
        
        MatrixPairwiseDist=nan(length(SylsUnique), length(SylsUnique));
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' exptname]);
        hold on;
        
        % === get pairwise between all syls
        for j=1:length(SylsUnique)
            syl1=SylsUnique{j};
            
            for jj=j+1:length(SylsUnique);
                
                syl2=SylsUnique{jj};
                
                disp([syl1 '-' syl2])
                
                % ==== collect acoustic distance
                fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl1).fv_baseline_PCAscore_mean;
                fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_PCAscore_mean;
                
                euclDist=sqrt(sum((fv1-fv2).^2));
                
                % ===== COLLECT OUTPUTS
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).ACOUSTIC.distanceFrom_using8dimPCA.(syl2)=euclDist;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).ACOUSTIC.distanceFrom_using8dimPCA.(syl1)=euclDist;
                
                AllPairwiseDist=[AllPairwiseDist euclDist];
                MatrixPairwiseDist(j,jj)=euclDist;
                
%                 lt_plot_text(j, jj, num2str(euclDist), 'b');
                
            end
        end
        
        % == plot heat map of pairwise distances
        imagesc(MatrixPairwiseDist, clims);
        colormap('hot');
        colorbar;
        
        CellsUnderThresh=MatrixPairwiseDist<=threshold;
        
        set(gca, 'XTick', 1:length(SylsUnique));
        set(gca, 'XTickLabel', SylsUnique);
        
        set(gca, 'YTick', 1:length(SylsUnique));
        set(gca, 'YTickLabel', SylsUnique);
        
        rotateXLabels(gca, 90);
        axis tight;
    end
end


%% == plot

lt_figure; 
lt_subplot(4,2,1); hold on;
lt_plot_histogram(AllPairwiseDist)

lt_subplot(4,2,2); hold on;
Y=AllPairwiseDist;
[f, xi, bw]=ksdensity(Y, 'bandwidth', bwidth);
plot(xi, f, '.-');
% plot line for threshold
line([threshold threshold], ylim);

lt_subplot(4,2,3); hold on;
xlabel('Acoustic Distance (z-scores');
ylabel('Probability');
Y=AllPairwiseDist;
[f, xi, bw]=ksdensity(Y, 'bandwidth', bwidth);
lt_plot(xi, f, {'Marker','none', 'LineStyle','-'});
line([threshold threshold], ylim, 'Color','k');





