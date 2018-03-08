%% lt - quick script for extracting specific bird and motifs into Motifstruct
% ====== by default does time warping.

close all; clear MOTIFSTATS_Compiled;
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 1;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;

% --- to make sure extracts motifs
MotifsToCollect = {'pu69wh78', {'(j)jjbhhg', '(a)abhhg'}};
MotifsToCollect = {'wh44wh39', {'(n)hh', '(d)kccbb'}};
Params_regexp.motif_predur = 0.05;
Params_regexp.motif_postdur = 0.15;
Params_regexp.preAndPostDurRelSameTimept = 0;
Params_regexp.RemoveIfTooLongGapDur = 1;

MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF, MotifsToCollect, Params_regexp);


if (0)
    %% ======================= LINEAR TIME WARP
    TimeWarpParams = {'pu69wh78', '(j)jjbhhg', [1:13], ...
        'pu69wh78', '(a)abhhg', [1:11]};
    NumBirds = length(MOTIFSTATS_Compiled.birds);
    
    MOTIFSTATS_Compiled.TimeWarpParams = TimeWarpParams;
    for i=1:NumBirds
        
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        
        nneur = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        for ii=1:nneur
            nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
            for iii=1:nummotifs
                
                motifthis = motiflist{iii};
                birdthis = MOTIFSTATS_Compiled.birds(i).birdname;
                
                ind1 = find(strcmp(TimeWarpParams, birdthis));
                ind2 = find(strcmp(TimeWarpParams, motifthis));
                
                ind3 = ind1(ind1 == ind2-1)+2; % actual ind of params
                
                if isempty(ind3)
                    % then this bird or motif is not specificed in params
                    disp(['PROBLEM - b ' birdname '- motif ' num2str(motifthis) ' NOT SPECIFIED (WILL NOT TIME WARP)']);
                    continue
                end
                
                disp([birdthis '-n' num2str(ii) '-' motifthis]);
                
                % ================= DO TIME WARP
                regionstowarp = TimeWarpParams{ind3};
                segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
                
                expectedsegs = 2*length(segextract(1).matchlabel) - 1;
                segextract = lt_neural_LinTimeWarpSegmented(segextract, ...
                    regionstowarp, expectedsegs);
                
                MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract = ...
                    segextract;
            end
        end
    end
    
else
    %% ==================== LINEAR TIME WARP, COMBINGIN ALL NEURONS
    TimeWarpParams = {'pu69wh78', '(j)jjbhhg', [1:13], ...
        'pu69wh78', '(a)abhhg', [1:11], ...
        'wh44wh39', '(n)hh', [1:5], ...
        'wh44wh39', '(d)kccbb', [1:11]};
    NumBirds = length(MOTIFSTATS_Compiled.birds);
    
    MOTIFSTATS_Compiled.TimeWarpParams = TimeWarpParams;
    
    for i=1:NumBirds
        
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        nummotifs = length(motiflist);
        nneur = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        
        for iii=1:nummotifs
            
            motifthis = motiflist{iii};
            birdthis = MOTIFSTATS_Compiled.birds(i).birdname;
            
            ind1 = find(strcmp(TimeWarpParams, birdthis));
            ind2 = find(strcmp(TimeWarpParams, motifthis));
            ind3 = ind1(ind1 == ind2-1)+2; % actual ind of params
            
            if isempty(ind3)
                % then this bird or motif is not specificed in params
                disp(['PROBLEM - b ' birdname '- motif ' num2str(motifthis) ' NOT SPECIFIED (WILL NOT TIME WARP)']);
                continue
            end
            
            %         for ii=1:nneur
            %             disp([birdthis '-n' num2str(ii) '-' motifthis]);
            %
            %             % ================= DO TIME WARP
            %             regionstowarp = TimeWarpParams{ind3};
            %             segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
            %
            %             expectedsegs = 2*length(segextract(1).matchlabel) - 1;
            %             segextract = lt_neural_LinTimeWarpSegmented(segextract, ...
            %                 regionstowarp, expectedsegs);
            %
            %             MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract = ...
            %                 segextract;
            %         end
            %
            % ================== combine segextract across all neurons
            segextractAll = [];
            segIndsAll = cell(1, nneur);  % --- save inds to be able to put back into specific neurons
            
            for ii=1:nneur
                disp([birdthis '-n' num2str(ii) '-' motifthis]);
                
                % =================
                segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
                segIndsAll{ii} = length(segextractAll)+1:(length(segextractAll)+length(segextract));
                
                % --- add to concated
                segextractAll = [segextractAll segextract];
            end
            
            % ============== DO TIMEWARP
            regionstowarp = TimeWarpParams{ind3};
            
            expectedsegs = 2*length(segextract(1).matchlabel) - 1;
            segextractAll = lt_neural_LinTimeWarpSegmented(segextractAll, ...
                regionstowarp, expectedsegs);
            
            % =========== SLIDE BACK INTO MOTIFSTATS
            for ii=1:nneur
                
                indstotake = segIndsAll{ii};
                assert(length(indstotake) == length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract), 'asfsda');
                
                MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract = ...
                    segextractAll(indstotake);
            end
            
            
        end
    end
end