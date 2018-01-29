function MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt, ...
    collectFF, MotifsToCollect, Params_regexp)
%% 1/17/18 - added ability to manually enter motifs you'd liek to collect
% ALSO ADDED ENTRY FOR REGEXP PARAMS (if no exist then uses defaults)


if ~exist('Params_regexp', 'var')
    Params_regexp = [];
end


if ~exist('MotifsToCollect', 'var')
    MotifsToCollect = [];
end

%%

if ~exist('collectFF', 'var')
    collectFF = 1;
end
if isempty(collectFF)
    collectFF=1;
end



%% lt 12/20/17 - allows you to organize by neurons, instead of by expt.
% NOTE: default (by expt) makes neuron numbers innacurate (i.e. they become
% numbered 1...N where N is number of neurons in this expt

if ~exist('OrganizeByExpt', 'var')
    OrganizeByExpt =1;
end
if isempty(OrganizeByExpt)
    OrganizeByExpt=1;
end
   


%% modified 10/16/17 by LT to save output struct. also asks before running
% if want to load old struct - finds structs that have same params as
% current params

if ~exist('saveOn', 'var')
    saveOn = 1;
end

if ~exist('onlyCollectTargSyl', 'var')
    onlyCollectTargSyl = 0;
end
if isempty(onlyCollectTargSyl)
    onlyCollectTargSyl=0;
end

%%
% LearnKeepOnlyBase for learning, if 1, then keeps only baseline periods (i.e. before WN)

%%

if ~exist('collectWNhit', 'var')
    collectWNhit = 1;
end

if ~exist('LearnKeepOnlyBase', 'var')
    LearnKeepOnlyBase =0;
end


%% ================== first check if this is done previously and saved
savedir = '/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled';
cd(savedir);
loadedold = 0;

% === check each params, keep list of those that are same as current params
prmslist = dir('Params_*');
for i=1:length(prmslist)
    
    load(prmslist(i).name);
    
    matchvec = [isempty(setxor(Params.BirdsToKeep, SummaryStruct.loadparams.BirdsToKeep)) ...
        isempty(setxor(Params.BrainArea, SummaryStruct.loadparams.BrainArea)) ...
        isempty(setxor(Params.ExptToKeep, SummaryStruct.loadparams.ExptToKeep)) ...
        isempty(setxor(Params.BatchesDesired, SummaryStruct.loadparams.BatchesDesired)) ...
        isempty(setxor(Params.RecordingDepth, SummaryStruct.loadparams.RecordingDepth)) ...
        isempty(setxor(Params.ChannelsDesired, SummaryStruct.loadparams.ChannelsDesired)) ...
        Params.LearningOnly == SummaryStruct.loadparams.LearningOnly ...
        Params.collectWNhit == collectWNhit ...
        Params.LearnKeepOnlyBase == LearnKeepOnlyBase];
    
    if all(matchvec==1)
        % --- this is a match
        disp(prmslist(i).name);
        disp(Params)
        if strcmp(input(['FOUND PREVIOUS EXACT MATCH - LOAD? extracted on ' num2str(Params.extractiontime) ' (y or n) '], 's'), 'y')
            
            
            
            indtmp = strfind(prmslist(i).name, '_');
            structname = ['MOTIFSTATS_Compiled_' prmslist(i).name(indtmp+1:end)];
            load(structname);
            loadedold = 1;
            disp('--- LOADED OLD MOTIFSTATS_Compiled !!');
            break
        elseif strcmp(input(['DELETE? (y or n) '], 's'), 'y')
                    indtmp = strfind(prmslist(i).name, '_');
            structname = ['MOTIFSTATS_Compiled_' prmslist(i).name(indtmp+1:end)];
    eval(['!rm ' structname]);
        end
        
    end
end



%% if for learning only keeping base, then expand the motifs to include all
% also get motifs if manually get motifs.

if ~isempty(MotifsToCollect) | LearnKeepOnlyBase==1
            SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct, LearnKeepOnlyBase, MotifsToCollect);
end



%%
if loadedold==0
    %% lt 6/8/17 - multiple birds, extracts motifs/segments, etc
    
    NumBirds = length(SummaryStruct.birds);
    MOTIFSTATS_Compiled = struct; % heirarchy: birds --> expt --> neurons
    
    % 
    
    for i=1:NumBirds
        
        birdname = SummaryStruct.birds(i).birdname;
        
        if OrganizeByExpt==1
            ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
            
            for ll=1:length(ListOfExpts)
                exptname = ListOfExpts{ll};
                
                inds = strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname);
                
                SummaryStruct_tmp = struct;
                SummaryStruct_tmp.birds(1).neurons = SummaryStruct.birds(i).neurons(inds);
                SummaryStruct_tmp.birds(1).birdname = SummaryStruct.birds(i).birdname;
                
                FFparams.collectFF=collectFF;
                
                % === extract for just this expt
                [MOTIFSTATS] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct_tmp, ...
                    collectWNhit, onlyCollectTargSyl, LearnKeepOnlyBase, FFparams, Params_regexp);
                
                % === OUTPUT
                MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS = MOTIFSTATS;
                MOTIFSTATS_Compiled.birds(i).exptnum(ll).neurIDOriginal_inorder = inds;
                MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct = SummaryStruct_tmp;
                MOTIFSTATS_Compiled.birds(i).exptnum(ll).exptname = exptname;
                MOTIFSTATS_Compiled.birds(i).birdname = birdname;
                
            end
        else
            
                SummaryStruct_tmp = struct;
                SummaryStruct_tmp.birds(1) = SummaryStruct.birds(i);
                
                FFparams.collectFF=collectFF;
                
                % === extract
                [MOTIFSTATS] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct_tmp, ...
                    collectWNhit, onlyCollectTargSyl, LearnKeepOnlyBase, FFparams, Params_regexp);
                
                % === OUTPUT
                MOTIFSTATS_Compiled.birds(i).MOTIFSTATS = MOTIFSTATS;
                MOTIFSTATS_Compiled.birds(i).SummaryStruct = SummaryStruct_tmp;
                MOTIFSTATS_Compiled.birds(i).birdname = birdname;
                MOTIFSTATS_Compiled.birds(i).Params_regexp = Params_regexp;

        end
    end
    
    %% --- save
    
    if saveOn == 1
        
        % ==== make a params structure
        Params = SummaryStruct.loadparams;
        Params.collectWNhit = collectWNhit;
        Params.LearnKeepOnlyBase = LearnKeepOnlyBase;
       Params.OrganizeByExpt = OrganizeByExpt;
        
        % ==== save
        cd(savedir)
        
        tstamp = lt_get_timestamp(0);
        
        save(['Params_' tstamp], 'Params');
        save(['MOTIFSTATS_Compiled_' tstamp], 'MOTIFSTATS_Compiled', '-v7.3');
        
        
    end
end
