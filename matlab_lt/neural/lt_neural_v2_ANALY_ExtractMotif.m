function [MOTIFSTATS, SummaryStruct] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct, ...
    collectWNhit, onlyCollectTargSyl, LearnKeepOnlyBase)
%% ONLY WORKS FOR SINGLE BIRD!!! - extracts motif information into one structure
% NOT LEARNING SPECIFIC, GENERAL USE

%% PARAMS

motif_predur = 0.15;
motif_postdur = 0.05;

% for segments extract:
alignByOnset=1;
WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
% those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
% +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error



if ~exist('collectWNhit', 'var')
    collectWNhit=1;
end

if ~exist('onlyCollectTargSyl', 'var')
    onlyCollectTargSyl=0; % if 1, then has to be learning experiment and so has targ syl defined.
end

if ~exist('LearnKeepOnlyBase', 'var');
    LearnKeepOnlyBase=0; % if 1, then keeps only baseline if you're a learning expt
end

%% what are motifs for this bird?
NumBirds = length(SummaryStruct.birds);
assert(NumBirds ==1, 'too many birds');

MotifsActual = SummaryStruct.birds(1).neurons(1).POSTINFO.MotifsActual;
motif_regexpr_str = SummaryStruct.birds(1).neurons(1).POSTINFO.MotifsActual_regexpStr;
singlesyls = SummaryStruct.birds(1).neurons(1).POSTINFO.SingleSyls_unique;


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


% motif_regexpr_str = {}; % will put motifs for extraction here
%
% for i=1:length(MotifsActual)
%     motif_actual = MotifsActual{i};
%
%     % for each vocalization in the motif, extract one segment
%     numvocals = length(motif_actual);
%
%     for ii=1:numvocals
%         if ii==1
%             segmentmotif = ['(' motif_actual(ii) ')' motif_actual(ii+1:end)];
%
%         elseif ii == numvocals
%             segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')'];
%
%         else
%
%             segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')' motif_actual(ii+1:end)];
%
%         end
%
%         %        disp(segmentmotif);
%         motif_regexpr_str = [motif_regexpr_str segmentmotif];
%     end
% end

%% extract data

% one struct for each neuron x syl
MOTIFSTATS = struct;

NumNeurons = length(SummaryStruct.birds(1).neurons);
NumSyls = length(motif_regexpr_str);
for i=1:NumNeurons
    cd(SummaryStruct.birds(1).neurons(i).dirname);
    
    exptname = SummaryStruct.birds(1).neurons(i).exptID;
    birdname =SummaryStruct.birds(1).birdname;
    % -- load data for this neuron
    batchf=SummaryStruct.birds(1).neurons(i).batchfilename;
    channel_board=SummaryStruct.birds(1).neurons(i).channel;
    if collectWNhit==0
        extractSound = 0;
    else
        extractSound = 1;
    end
    cd ..
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractSound);
            Params.birdname = SummaryStruct.birds(1).birdname;
        Params.exptname = SummaryStruct.birds(1).neurons(i).exptID;

    % === if only extract targ syl, find what targ syl is
    
    if onlyCollectTargSyl==1
        
        tmp = lt_neural_v2_LoadLearnMetadat;
        indbird = strcmp({tmp.bird.birdname}, birdname);
        indexpt = strcmp(tmp.bird(indbird).info(1,:), exptname);
        TargSyls = tmp.bird(indbird).info(2,indexpt);
        assert(~isempty(TargSyls), 'PROBELM< DEFINE TARG SYLS FIRST'); % must not be empty
    end
    
    
    % --- extract data, once for each motif
    for j=1:NumSyls
        regexpr_str=motif_regexpr_str{j};
        
        if onlyCollectTargSyl==1
            % make sure this is a targ syl
            if any(strcmp(regexpr_str, TargSyls))
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                    regexpr_str, motif_predur, motif_postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
                    0, 1, collectWNhit, 1, LearnKeepOnlyBase);
                
            else
                SegmentsExtract=struct;
                Params=struct;
            end
        else
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                regexpr_str, motif_predur, motif_postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
                0, 1, collectWNhit, 1, LearnKeepOnlyBase);
            
        end
        
        %         [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
        %             regexpr_str, motif_predur, motif_postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
        %             0, 1, collectWNhit);
        %         close all;
        
        
        %%  -- throw out all data post switchtime [obsolete - doing in regexp code now]
%     [islearning, LearnSummary, switchtime] = lt_neural_v2_QUICK_islearning(birdname, exptname, 1);
%         if isfield(SegmentsExtract, 'song_datenum')
%             % i.e., has data
%             
%         if LearnKeepOnlyBase ==1
%             if islearning==1
%                 assert(~isempty(switchtime), 'asdfasdf')
%                 
%                 tvals = [SegmentsExtract.song_datenum];
%                 indstoremove = tvals > switchtime;
%                 
%                 SegmentsExtract(indstoremove) = [];
%                 
%                 disp([' === removed WN files from (learning) : ' birdname '-' exptname '-neur' num2str(i) '-' regexpr_str ' - ' num2str(sum(indstoremove)) '/' num2str(length(indstoremove))]);
%                 
%                 
%                 % --- if all data removed, make it an empty structure
%                 if all(indstoremove)
%                     SegmentsExtract = struct;
%                 end
%             end
%         end
%         end
        
        
        
        %% OUTPUT
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


%% ==== if this is learning neuron, choose to only keep baseline



%% other things to output
MOTIFSTATS.params.motif_predur = motif_predur;
MOTIFSTATS.params.motif_postdur = motif_postdur;
MOTIFSTATS.params.motif_regexpr_str = motif_regexpr_str;
MOTIFSTATS.params.singlesyls_unique = singlesyls;
% MOTIFSTATS.params.BirdsToKeep = BirdsToKeep;
% MOTIFSTATS.params.BrainArea = BrainArea;
% MOTIFSTATS.params.ExptToKeep = ExptToKeep;
MOTIFSTATS.params.MotifsActual = MotifsActual;
% MOTIFSTATS.params.TargSyls = TargSyls;




