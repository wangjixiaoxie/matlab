function MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn)
%% modified 10/16/17 by LT to save output struct. also asks before running
% if want to load old struct - finds structs that have same params as
% current params

if ~exist('saveOn', 'var')
    saveOn = 1;
end


%%
% LearnKeepOnlyBase for learning, if 1, then keeps only baseline periods (i.e. before WN)

%%

if ~exist('collectWNhit', 'var')
    collectWNhit = 1;
end

if ~exist('LearnKeepOnlyBase', 'var');
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
    
    matchvec = [isempty(setdiff(Params.BirdsToKeep, SummaryStruct.loadparams.BirdsToKeep)) ...
        isempty(setdiff(Params.BrainArea, SummaryStruct.loadparams.BrainArea)) ...
        isempty(setdiff(Params.ExptToKeep, SummaryStruct.loadparams.ExptToKeep)) ...
        isempty(setdiff(Params.BatchesDesired, SummaryStruct.loadparams.BatchesDesired)) ...
        isempty(setdiff(Params.RecordingDepth, SummaryStruct.loadparams.RecordingDepth)) ...
        isempty(setdiff(Params.ChannelsDesired, SummaryStruct.loadparams.ChannelsDesired)) ...
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
        end
        
    end
end


if loadedold==0
    %% lt 6/8/17 - multiple birds, extracts motifs/segments, etc
    
    NumBirds = length(SummaryStruct.birds);
    MOTIFSTATS_Compiled = struct; % heirarchy: birds --> expt --> neurons
    
    for i=1:NumBirds
        
        ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
        birdname = SummaryStruct.birds(i).birdname;
        
        for ll=1:length(ListOfExpts)
            exptname = ListOfExpts{ll};
            
            inds = strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname);
            
            SummaryStruct_tmp = struct;
            SummaryStruct_tmp.birds(1).neurons = SummaryStruct.birds(i).neurons(inds);
            SummaryStruct_tmp.birds(1).birdname = SummaryStruct.birds(i).birdname;
            
            % === extract for just this expt
            [MOTIFSTATS] = lt_neural_v2_ANALY_ExtractMotif(SummaryStruct_tmp, ...
                collectWNhit, 0, LearnKeepOnlyBase);
            
            % === OUTPUT
            MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS = MOTIFSTATS;
            MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct = SummaryStruct_tmp;
            MOTIFSTATS_Compiled.birds(i).exptnum(ll).exptname = exptname;
            MOTIFSTATS_Compiled.birds(i).birdname = birdname;
            
        end
    end
    
    %% --- save
    
    if saveOn == 1;
        
        % ==== make a params structure
        Params = SummaryStruct.loadparams;
        Params.collectWNhit = collectWNhit;
        Params.LearnKeepOnlyBase = LearnKeepOnlyBase;
        
        
        % ==== save
        cd(savedir)
        
        tstamp = lt_get_timestamp(0);
        
        save(['Params_' tstamp], 'Params');
        save(['MOTIFSTATS_Compiled_' tstamp], 'MOTIFSTATS_Compiled', '-v7.3');
        
        
    end
end
