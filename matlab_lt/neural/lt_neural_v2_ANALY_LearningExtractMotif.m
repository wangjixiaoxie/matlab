function [MOTIFSTATS, SummaryStruct] = lt_neural_v2_ANALY_LearningExtractMotif(SummaryStruct)
%% PARAMS

motif_predur = 0.15;
motif_postdur = 0.05;


%% first extract neurons [actually do this before run this code]

% % === 1) same learning experiment
% % ===== LMANlearn2
% BirdsToKeep = {'wh6pk36'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {};
% ExptToKeep = {'LMANlearn2'};
% [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, BrainArea, ExptToKeep);


%% what are motifs for this bird?
NumBirds = length(SummaryStruct.birds);
assert(NumBirds ==1, 'too many birds');

MotifsActual = SummaryStruct.birds(1).neurons(1).POSTINFO.MotifsActual;

% === GIVEN ACTUAL MOTIFS FOR THIS BIRD, DETERMINE WHICH SYLS TO LOOK AT
if strcmp(SummaryStruct.birds(1).birdname, 'wh6pk36')
%     MotifsActual = {'nlcchb', 'jklcchb', 'ga', 'mksd' ,'vb'};
    TargSyls = {'nlcch(b)', 'jklcch(b)'};
elseif strcmp(SummaryStruct.birds(1).birdname, 'bk7')
%     MotifsActual = {'nnhh', 'gh', 'vbbb', 'gb', 'jkk', 'kl', 'gv', 'yoo', 'rs'}; % COULD IMPROVE
    TargSyls = {'g(h)'};
elseif strcmp(SummaryStruct.birds(1).birdname, 'bu77wh13')
%     MotifsActual = {'kspj', 'ab', 'bh', 'jb', 'nkrs'}; % COULD IMPROVE
    %       MotifsActual = {'kspj', 'ab', 'bh', 'ijbh', 'nkrs'}; % COULD IMPROVE
    TargSyls = {'a(b)'};
elseif strcmp(SummaryStruct.birds(1).birdname, 'br92br54')
%     MotifsActual = {'nkh', 'ddd', 'dh', 'agc', 'cc'}; % COULD IMPROVE
    %       MotifsActual = {'kspj', 'ab', 'bh', 'ijbh', 'nkrs'}; % COULD IMPROVE
    if strcmp(SummaryStruct.birds(1).neurons(1).exptID, 'LMANlearn2')
        TargSyls = {'nk(h)'};
    elseif strcmp(SummaryStruct.birds(1).neurons(1).exptID, 'LMANlearn3')
        TargSyls = {'ag(c)'};
    elseif strcmp(SummaryStruct.birds(1).neurons(1).exptID, 'LMANlearn4')
         TargSyls = {'nk(h)'};       
    elseif strcmp(SummaryStruct.birds(1).neurons(1).exptID, 'LMANlearn5')
         TargSyls = {'nk(h)', 'd(h)'};       
    end
elseif strcmp(SummaryStruct.birds(1).birdname, 'or74bk35')
         TargSyls = {'an(b)'};       
end


motif_regexpr_str = {}; % will put motifs for extraction here

for i=1:length(MotifsActual)
    motif_actual = MotifsActual{i};
    
    % for each vocalization in the motif, extract one segment
    numvocals = length(motif_actual);
    
    for ii=1:numvocals
        if ii==1
            segmentmotif = ['(' motif_actual(ii) ')' motif_actual(ii+1:end)];
            
        elseif ii == numvocals
            segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')'];
            
        else
            
            segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')' motif_actual(ii+1:end)];
            
        end
        
        %        disp(segmentmotif);
        
        motif_regexpr_str = [motif_regexpr_str segmentmotif];
        
    end
end



%% extract data

% one struct for each neuron x syl
MOTIFSTATS = struct;

NumNeurons = length(SummaryStruct.birds(1).neurons);
NumSyls = length(motif_regexpr_str);

for i=1:NumNeurons
    cd(SummaryStruct.birds(1).neurons(i).dirname);
    
    % -- load data for this neuron
    batchf=SummaryStruct.birds(1).neurons(i).batchfilename;
    channel_board=SummaryStruct.birds(1).neurons(i).channel;
    extractSound = 1;
    cd ..
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractSound);
    
    % --- extract data, once for each motif
    for j=1:NumSyls
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
        FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
        % +1 is 1 after token
        FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams);
        close all;
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
        
    end
    MOTIFSTATS.neurons(i).clustnum = SummaryStruct.birds(1).neurons(i).clustnum;
    
end


%% === MAKE SURE ALL TRIALS ARE IN TEMPORAL ORDER

for i=1:NumNeurons
    for m=1:NumSyls
        
        segextract = MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        
        if ~isfield(segextract, 'song_datenum')
        continue
        end
        
        all_datenums=[segextract.song_datenum];
        
        [~, inds] = sort(all_datenums);
        
        if all(diff(inds))~=1
            disp('HAD TO REORDER !! ---- not a problem');
            
            MOTIFSTATS.neurons(i).motif(m).SegmentsExtract=segextract(inds);
        end
    end
end

%% other things to output
MOTIFSTATS.params.motif_predur = motif_predur;
MOTIFSTATS.params.motif_postdur = motif_postdur;
MOTIFSTATS.params.motif_regexpr_str = motif_regexpr_str;
% MOTIFSTATS.params.BirdsToKeep = BirdsToKeep;
% MOTIFSTATS.params.BrainArea = BrainArea;
% MOTIFSTATS.params.ExptToKeep = ExptToKeep;
MOTIFSTATS.params.MotifsActual = MotifsActual;
MOTIFSTATS.params.TargSyls = TargSyls;




