close all;

%% EXTRACT FF AND SAVE

FFparamsAll.bird(1).birdname = 'bk7';
FFparamsAll.bird(1).FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
            'v', [2450 4300]};
FFparamsAll.bird(1).FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
FFparamsAll.bird(1).FFparams.cell_of_FFtimebins_DurLearn={'h', [0.034 0.038], 'b', [0.053 0.07], ...
            'v', [0.052 0.07]}; % WN on g H

       
FFparamsAll.bird(2).birdname = 'bu77wh13';
FFparamsAll.bird(2).FFparams.cell_of_freqwinds={'b', [2700 3900], 'h', [2600 3900], 'a', [1300 2600]};
FFparamsAll.bird(2).FFparams.cell_of_FFtimebins={'b', [0.03 0.038], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)
FFparamsAll.bird(2).FFparams.cell_of_FFtimebins_DurLearn={'b', [0.03 0.038], 'h', [0.038 0.052], 'a', [0.056 0.081]}; % in sec, relative to onset (i.e. see vector T)


FFparamsAll.bird(3).birdname = 'wh6pk36';
FFparamsAll.bird(3).FFparams.cell_of_freqwinds={'c', [2100 3100], 'h', [2800 4000], 'b', [2700 3800], ...
    'a', [1300 2200], 's', [4000 5100], 'd', [900 2000],  'n', [3300 4300], 'v', [2600 4000]};
FFparamsAll.bird(3).FFparams.cell_of_FFtimebins={'c', [0.035 0.05], 'h', [0.023 0.033], 'b', [0.034 0.038], ...
    'a', [0.06 0.08], 's', [0.02 0.023], 'd', [0.02 0.027],  'n', [0.029 0.05], 'v', [0.028 0.036]};
FFparamsAll.bird(3).FFparams.cell_of_FFtimebins_DurLearn={'c', [0.035 0.05], 'h', [0.023 0.033], 'b', [0.034 0.038], ...
    'a', [0.06 0.08], 's', [0.02 0.023], 'd', [0.02 0.027],  'n', [0.029 0.05], 'v', [0.028 0.036]};


% TEMPORARY
FFparamsAll.bird(4).birdname = 'br92br54';
FFparamsAll.bird(4).FFparams.cell_of_freqwinds={'a', [750 1400], 'c', [1200 1800], ...
    'h', [2350 3900], 'd', [1300 3400], 'k', [800 1800]};
FFparamsAll.bird(4).FFparams.cell_of_FFtimebins={'a', [0.057 0.08], 'c', [0.044 0.055], ...
    'h', [0.040 0.049], 'd', [0.032 0.052], 'k', [0.05 0.055]};
FFparamsAll.bird(4).FFparams.cell_of_FFtimebins_DurLearn={'a', [0.057 0.08], 'c', [0.044 0.055], ...
    'h', [0.040 0.049], 'd', [0.029 0.033], 'k', [0.05 0.055]};


FFparamsAll.bird(5).birdname = 'or74bk35';
FFparamsAll.bird(5).FFparams.cell_of_freqwinds={'a', [1100 2700], 'g', [1000 2300], ...
    'n', [3000 4300], 'b', [2750 3900]};
FFparamsAll.bird(5).FFparams.cell_of_FFtimebins={'a', [0.04 0.06], 'g', [0.095 0.105], ...
    'n', [0.04 0.05], 'b', [0.033 0.041]};
FFparamsAll.bird(5).FFparams.cell_of_FFtimebins_DurLearn={'a', [0.04 0.06], 'g', [0.095 0.105], ...
    'n', [0.04 0.05], 'b', [0.033 0.041]};

overWrite = 0; % note, will overwrite rgardless if detects chagnes (NOTE: always overwrites if detects changes)
plotSpec = 0; % to plot raw spec overlayed with PC and windows.
plotOnSong = 41; % will only start plotting spec once hit this song num.
plotSyl = ''; % to focus on just one syl. NOT DONE YET
equalizeParams = 0; % if 1, then makes sure t and f windows match. if 0 then only makes sure labels match.
lt_neural_v2_EXTRACT_FF_tmp(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl, equalizeParams);


%% =================== LEARNING 
% CAN DO MULTIPLE BIRDS
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit);


% =========== PICK OUT SAME TYPE/DIFF [LEANRING SPECIFIC]
numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:numbirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
   
    for ii=1:numexpts
       exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        % ---- EXTRACT TARG SYL
        tmp = lt_neural_v2_LoadLearnMetadat;
        indbird = strcmp({tmp.bird.birdname}, birdname);
        indexpt = strcmp(tmp.bird(indbird).info(1,:), exptname);
        TargSyls = tmp.bird(indbird).info(2,indexpt);

        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls = TargSyls;
        
        % ----
%         TargSyls = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls;
        MotifsActual = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.MotifsActual;
        
        [SameTypeSyls, DiffTypeSyls, motif_regexpr_str] = ...
        lt_neural_v2_extractSameType(MotifsActual, TargSyls);

    % --- OUTPUT
MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.SameTypeSyls = SameTypeSyls;
MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.DiffTypeSyls = DiffTypeSyls;

    end
end


%% ==== FOR EACH LEARNING EXPERIMENT, PLOT TIMELINE OF NEURONS AND LEARNING
NumBirds = length(SummaryStruct.birds);
for i=1:NumBirds
   ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
   
   for ll = 1:length(ListOfExpts)
        
       MOTIFSTATS = MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS;
       SummaryStruct_tmp = MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct;
       
       birdname = SummaryStruct.birds(i).birdname;
       exptname = ListOfExpts{ll};
       
       if ~isfield(MOTIFSTATS.params, 'TargSyls')
           disp('asdasdf');
       end
       if isempty(MOTIFSTATS.params.TargSyls);
           sdafasdf
       end
       
      % === PLOT
      lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_tmp, MOTIFSTATS);
       title([birdname '-' exptname]);
   end
end
   
% ========= SAVE AND CLOSE
lt_save_figs_to_folder('062117_AllLearn_Summary', 1);

