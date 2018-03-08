function [FRmatMotifByNeur, AllMotifRegexp, AllNeurLocation] = ...
    lt_neural_v2_PseudoPop_Master(MOTIFSTATS_Compiled, NormToNeurMean, ...
    premotorWind, plotnans)

%% ============= PARAMS

% NormToNeurMean = 1;
% premotorWind = [-0.025 0.035];
% % premotorWind = [0.01 0.07];
% % premotorWind = [-0.5 0.01];


%% =================== EXTRACT DATA

i=1;
numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);

motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
motif_predur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;

FRmatMotifByNeur = nan(length(motiflist), numneurons);


AllNeurLocation = {};
AllMotifRegexp = motiflist;

for ii=1:numneurons
    
    
    nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
    
    loc = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).NOTE_Location;
    AllNeurLocation = [AllNeurLocation; loc];
    
    for iii=1:nummotifs
        
        motifthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str{iii};
        assert(strcmp(motifthis, motiflist{iii})==1, 'fsdaf');
        
        segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
        
        
        % ====================== EXTRACT MEAN FR IN PREMOTOR WINDOW
        Nspks = lt_neural_QUICK_SpkCnts(segextract, motif_predur, ...
            premotorWind);
        Nspks_mean = mean(Nspks);
        Nspks_trialCV = std(Nspks)/mean(Nspks);
        Nspks_trialSD = std(Nspks);
        
        % ---------- divide into 2 halves to get higher dimension
        % representation.
        
        % ====================== PUT INTO CORRECT SLOT (MOTIF, NEURON)
        FRmatMotifByNeur(iii, ii) = Nspks_mean;
    end
end

%% ====================== pare down matrix is any missing values.
if any(isnan(FRmatMotifByNeur(:)))
    if plotnans==1
    lt_figure; hold on;
    
    % ------------ THEN WILL HAVE TO PARE DOWN MATRIX
    % ---- 1) plot matrix before paring down
    lt_subplot(2,2,1); hold on;
    xlabel('neurn');
    ylabel('motif');
    title('before paring down matrix (dot = empty)');
    spy(~isnan(FRmatMotifByNeur));
    end
    
    % ----- 2) remove neurons that are bad
    nummotifs = size(FRmatMotifByNeur,1);
    colsToRemove = sum(isnan(FRmatMotifByNeur),1) > nummotifs/4;
    
    FRmatMotifByNeur(:,colsToRemove) = [];
    AllNeurLocation(colsToRemove) = [];
    
    if plotnans==1
    lt_subplot(2,2,2); hold on;
    title('step1: after removing neurons missing 1/4 motifs');
    spy(~isnan(FRmatMotifByNeur));
    end
    
    % ---------3 now remove motifs that are missing neurons
    rowsToRemove = any(isnan(FRmatMotifByNeur)');
    FRmatMotifByNeur(rowsToRemove, :) = [];
    AllMotifRegexp(rowsToRemove) = [];
    
    if plotnans==1
    lt_subplot(2,2,3); hold on;
    title('step2: after removing motifs that miss any neurons');
    spy(~isnan(FRmatMotifByNeur));
    end
end


%% ===================== FOR EACH NEURON, SUBTRACT ITS OWN MEAN(ACROSS
% MOTIFS)
if NormToNeurMean==1
    
    tmp = repmat(mean(FRmatMotifByNeur,1), size(FRmatMotifByNeur,1), 1);
    FRmatMotifByNeur = FRmatMotifByNeur - tmp;
end




